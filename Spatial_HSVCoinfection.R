# Spatial model of HSV coinfection in skin (2D) - List-based version
# Pavitra Roychoudhury
# Fred Hutchinson Cancer Research Center
# Tracks two virus types and coinfection without recombination

# ini ---------------------------------------------------------------------

rm(list = ls()); 
library(RColorBrewer);

# GLOBAL SETTINGS ---------------------------------------------------------

SETTINGS <- list()

SETTINGS$L              <- 200;             # Specify grid dimensions
SETTINGS$dirname        <- 'Results/';   	# Directory to store results
SETTINGS$out_fname      <- 'results.tsv';	# File to store cell and virus counts
SETTINGS$sim_length     <- 1;				# Length of simulation in hours
SETTINGS$save_grids     <- TRUE;            # Save grids at the end of the sim? 
SETTINGS$visualize      <-'off';            # Options: 'screen','pdf','off','gif'
SETTINGS$refr_freq      <-25;				# Refresh frequency (hrs)
SETTINGS$site_scale     <-50;               # how many microns is one site
SETTINGS$place_virus1   <-list(type='infected producer',num=10);  #Initial virus 1
SETTINGS$place_virus2   <-list(type='infected producer',num=10);  #Initial virus 2

# GLOBAL PARAMETERS -------------------------------------------------------

PARAMETERS<-list(); 

# Rates (per site or cell per day)
PARAMETERS$max_rate                 <-25;       # Sets time resolution, should be greater than or equal to sum of all rates for a given state

# Uninfected cells
PARAMETERS$uninf_cell_death         <-0;        # Rate of uninfected cell death 
PARAMETERS$infectivity_free         <-0.1;      # Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
PARAMETERS$cell_div                 <-0.077;    # Rate of generation of new susceptible cells (repop of empty sites)

# Infected cells
PARAMETERS$inf_cell_death           <-1.25;     # Rate of infected cell death due to infection
PARAMETERS$viral_lag                <-8;		# Rate at which inf non-producer becomes viral producing cell (per cell per day)
PARAMETERS$vir_prod_rate_mean       <-100000;	# Rate at which an infected cell produces free virus (per cell per day, mean and sd)
PARAMETERS$vir_prod_rate_sd         <-10000;

# Virus
PARAMETERS$diff_rate_virus          <- 12;      # Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
PARAMETERS$freevirus_decay          <- 8.8;     # Rate of decay of free virus (per day)
PARAMETERS$cavirus_decay            <- 8.8;     # Rate of decay of cell-associated virus (per day)
PARAMETERS$viral_ejection_rate      <- 10;      # Rate of conversion from cell-associated to free virus (per day)
PARAMETERS$ca_to_free_proportion    <- 0.9;     # Proportion of CA virus converted to free virus per ejection event

# Probabilities, fractions, etc
PARAMETERS$host_density             <-1;        # Fraction of sites occupied by susceptible cells at start	
PARAMETERS$diff_range               <-3;        # n x n site, default n = 3, i.e. immediate Moore neighbourhood. 

SETTINGS$tot_updates <- SETTINGS$sim_length*PARAMETERS$max_rate; #Compute simulation length in number of updates


# setup -------------------------------------------------------------------

setup <- function() {
  #Directories
  print(paste('Results directory =',SETTINGS$dirname));
  dir.create(SETTINGS$dirname,showWarnings=FALSE);
  save(PARAMETERS,file=paste(SETTINGS$dirname,'/simulationparameters',sep=''));
  save(SETTINGS,file=paste(SETTINGS$dirname,'/simulationsettings',sep=''));
  if(SETTINGS$visualize=='gif') dir.create(paste(SETTINGS$dirname,'/imgs',sep=''),showWarnings=FALSE)
  
  #Dynamic variables - this stores time elapsed and population counts
  dynamic <<- list(num_updates=0,time=0,
                   susceptible_cells=0,infected_nonprodcells=0,infected_prodcells=0,dead_cells=0,
                   virions1=0,ca_virions1=0,virions2=0,ca_virions2=0,
                   max_vl=0,max_vl_time=0,
                   viral_cells1=0,ca_viral_cells1=0,viral_cells2=0,ca_viral_cells2=0,
                   coinfected_cells=0);
  output_results(dynamic,T); #write header
  
  #Visualization setup
  if(SETTINGS$visualize!='off'){
      if(SETTINGS$visualize=='pdf'){
          pdf(file=paste(SETTINGS$dirname,'/Results.pdf',sep=''),
              width=8,height=4); 
      }else if(SETTINGS$visualize=='screen'){
          quartz(width=8,height=4);   
      }else if(SETTINGS$visualize=='gif'){
          # library(animation)
          # png(file=paste(settings$dirname,'/imgs/Results.png',sep=''))
      }else{
          stop('Not a valid setting for "visualize"')
      }   	
      par(mar=c(1,1,2,1),oma=c(1,0,0,0));
      layout(matrix(c(1:2),1,2,byrow = TRUE),widths=c(1,1),heights=c(1,1));
  }
  gc(); 
}

# main --------------------------------------------------------------------
Spatial_main <- function() {
    print('Starting simulation.. hold on to your hat!');
    setup();
    init_grid();
    simulation();
}

# helper functions --------------------------------------------------------

# 2. Initialize grid and other arrays
init_grid <- function() {

    #Create spatial grid as an environment for fast hashed access
    #Each grid location is accessed via key "i,j" and contains a list with 
    # cell_state and virus counts for each type
    
    # spatial_grid <<- new.env(hash=TRUE);
    # 
    # for(i in 1:SETTINGS$L) {
    #     for(j in 1:SETTINGS$L) {
    #         key <- paste(i,j,sep=",");
    #         
    #         spatial_grid[[key]] <<- list(
    #             cell_state=as.integer(0),
    #             free_virus1=as.integer(0), 
    #             ca_virus1=as.integer(0),
    #             free_virus2=as.integer(0),
    #             =as.integer(0));
    #     }
    # }
    
    spatial_grid <<- list()
    spatial_grid$cell_state <<- matrix(0L, nrow = SETTINGS$L, ncol = SETTINGS$L)
    spatial_grid$free_virus1 <<- matrix(0L, nrow = SETTINGS$L, ncol = SETTINGS$L)
    spatial_grid$ca_virus1 <<- matrix(0L, nrow = SETTINGS$L, ncol = SETTINGS$L)
    spatial_grid$free_virus2 <<- matrix(0L, nrow = SETTINGS$L, ncol = SETTINGS$L)
    spatial_grid$ca_virus2 <<- matrix(0L, nrow = SETTINGS$L, ncol = SETTINGS$L)
    
}

# 3. SIMULATION
simulation<-function(){
  
  gaussian_mask_diffusion<-compute_gaussian(PARAMETERS$diff_range); 
  
  place_host(); #Place host cells
  if(SETTINGS$place_virus1$type!='none') place_virus(virus_type=1); #Place virus 1
  if(SETTINGS$place_virus2$type!='none') place_virus(virus_type=2); #Place virus 2
  
  dynamic$num_updates<<-0;
  compute_totals();
  output_results(dynamic,F);
  if(SETTINGS$visualize!='off'){
    visualize(); 
  }
  next_screen_grab<-dynamic$num_updates+SETTINGS$refr_freq;
  
  while(dynamic$num_updates < SETTINGS$tot_updates){
    
    #Generate random order of sites to update
    update_order<-sample(seq(1,SETTINGS$L^2),SETTINGS$L^2,replace=TRUE)
    row_order<-((update_order-1)%%SETTINGS$L)+1;
    col_order<-floor((update_order-1)/SETTINGS$L)+1;
    dice_cells_all<-runif(SETTINGS$L^2,min=0,max=1);
    dice_v1_all<-runif(SETTINGS$L^2,min=0,max=1);
    dice_v2_all<-runif(SETTINGS$L^2,min=0,max=1);
    dice_cav1_all<-runif(SETTINGS$L^2,min=0,max=1);
    dice_cav2_all<-runif(SETTINGS$L^2,min=0,max=1);
    
    #Update sites in random order (same as picking a random site to update one at a time)
    for(rep in 1:SETTINGS$L^2){
      i<-row_order[rep];
      j<-col_order[rep];
      dice_cells<-dice_cells_all[rep];
      dice_v1<-dice_v1_all[rep];
      dice_v2<-dice_v2_all[rep];
      dice_cav1<-dice_cav1_all[rep];
      dice_cav2<-dice_cav2_all[rep];
      
      #key<-paste(i,j,sep=",");
      #site  <- spatial_grid[[key]];
      site <- list(
          cell_state = spatial_grid$cell_state[i, j],
          free_virus1 = spatial_grid$free_virus1[i, j],
          ca_virus1 = spatial_grid$ca_virus1[i, j],
          free_virus2 = spatial_grid$free_virus2[i, j],
          ca_virus2 = spatial_grid$ca_virus2[i, j]
      )
      
      ##Host cells
      ## 0: vacant site
      ## 1: susceptible keratinocyte
      ## 2: infected, not producing virus
      ## 3: infected, producing virus
      ## 4: dead keratinocyte
      ##--------------------------------------------
      #State 0 or 4: vacant site or occupied by dead cell
      if(site$cell_state==0||site$cell_state==4){
        
        #Event 0.1: occupation by susceptible
        if(dice_cells<=PARAMETERS$cell_div/PARAMETERS$max_rate){
          #spatial_grid[[key]]$cell_state<-1;
          spatial_grid$cell_state[i, j] <<- 1;
        }
        
        
        #State 1: uninfected cell 
      }else if(site$cell_state==1){
        
        #Event 1.1: infection due to free virus at the site (any virus type can infect)
        total_free_virus<-site$free_virus1+site$free_virus2;
        if(total_free_virus>0 && dice_cells<=PARAMETERS$infectivity_free*total_free_virus/PARAMETERS$max_rate){
          #spatial_grid[[key]]$cell_state<-2; #becomes an infected cell
          spatial_grid$cell_state[i, j] <<- 2;
          #Remove one virus of random type
          virus_probs<-c(site$free_virus1,site$free_virus2)/total_free_virus;
          infecting_type<-sample(1:2,1,prob=virus_probs);
          if(infecting_type==1){
            #spatial_grid[[key]]$free_virus1<-site$free_virus1-1;
            spatial_grid$free_virus1[i, j] <<- site$free_virus1 - 1
          }else{
            #spatial_grid[[key]]$free_virus2<-site$free_virus2-1;
            spatial_grid$free_virus2[i, j] <<- site$free_virus2 - 1
          }
          
          #Event 1.2: apoptosis
        }else if(dice_cells<=PARAMETERS$uninf_cell_death/PARAMETERS$max_rate){
          #spatial_grid[[key]]$cell_state<-4;
          spatial_grid$cell_state[i, j] <<- 4
        } 
        
        #State 2: infected cell (non-virus producing)
      }else if(site$cell_state==2){					
        
        #Event 2.1: become a virus-producing cell
        if(dice_cells<=PARAMETERS$viral_lag/PARAMETERS$max_rate){
          #spatial_grid[[key]]$cell_state<-3; #becomes virus-producing infected cell
          spatial_grid$cell_state[i, j] <<- 3
          
          #Event 2.2: dies 
        }else if(dice_cells<=(PARAMETERS$viral_lag+PARAMETERS$uninf_cell_death)/PARAMETERS$max_rate){
          #spatial_grid[[key]]$cell_state<-4; #death of infected cell
          spatial_grid$cell_state[i, j] <<- 4
        }
        
        #State 3: infected cell (producing virus)
      }else if(site$cell_state==3){
        
        #Event 3.1: die
        if(dice_cells<=PARAMETERS$inf_cell_death/PARAMETERS$max_rate){
          #spatial_grid[[key]]$cell_state<-4; #death of infected cell
          spatial_grid$cell_state[i, j] <<- 4
          
          #Event 3.2: produce cell-associated virus (virus 1)
        }else{
          new_virus1<-max(0,round(rnorm(1,mean=PARAMETERS$vir_prod_rate_mean/PARAMETERS$max_rate,
                                        sd=PARAMETERS$vir_prod_rate_sd/PARAMETERS$max_rate)));
          #spatial_grid[[key]]$ca_virus1<-site$ca_virus1+new_virus1;
          spatial_grid$ca_virus1[i, j] <<- site$ca_virus1 + new_virus1
          
          #Event 3.3: produce cell-associated virus (virus 2)
          new_virus2<-max(0,round(rnorm(1,mean=PARAMETERS$vir_prod_rate_mean/PARAMETERS$max_rate,
                                        sd=PARAMETERS$vir_prod_rate_sd/PARAMETERS$max_rate)));
          #spatial_grid[[key]]$ca_virus2<-site$ca_virus2+new_virus2;
          spatial_grid$ca_virus2[i, j] <<- site$ca_virus2 + new_virus2
        }
      }
      
      
      ##Cell-associated virus - virus 1
      ##--------------------------------------------
      #site<-spatial_grid[[key]];  #Re-fetch site after cell state updates
      site <- list(
          cell_state = spatial_grid$cell_state[i, j],
          free_virus1 = spatial_grid$free_virus1[i, j],
          ca_virus1 = spatial_grid$ca_virus1[i, j],
          free_virus2 = spatial_grid$free_virus2[i, j],
          ca_virus2 = spatial_grid$ca_virus2[i, j]
      )
      
      if(site$ca_virus1>0){
        
        #Event 0a: Decay
        ca_viral_decay(i,j,dice_cav1,virus_type=1);
        
        #Event 0b: Ejection to free virus (probability increases with number of CA viruses)
        #site<-spatial_grid[[key]];  #Re-fetch after decay
        site <- list(
            cell_state = spatial_grid$cell_state[i, j],
            free_virus1 = spatial_grid$free_virus1[i, j],
            ca_virus1 = spatial_grid$ca_virus1[i, j],
            free_virus2 = spatial_grid$free_virus2[i, j],
            ca_virus2 = spatial_grid$ca_virus2[i, j]
        )
        if(site$ca_virus1>0 && dice_cav1<=PARAMETERS$viral_ejection_rate*site$ca_virus1/PARAMETERS$max_rate){
          ca_viral_ejection(i,j,virus_type=1);
        }
      }
      
      ##Cell-associated virus - virus 2
      ##--------------------------------------------
      #site<-spatial_grid[[key]];
      site <- list(
          cell_state = spatial_grid$cell_state[i, j],
          free_virus1 = spatial_grid$free_virus1[i, j],
          ca_virus1 = spatial_grid$ca_virus1[i, j],
          free_virus2 = spatial_grid$free_virus2[i, j],
          ca_virus2 = spatial_grid$ca_virus2[i, j]
      )
      if(site$ca_virus2>0){
        
        #Event 0a: Decay
        ca_viral_decay(i,j,dice_cav2,virus_type=2);
        
        #Event 0b: Ejection to free virus
        #site<-spatial_grid[[key]];
        site <- list(
            cell_state = spatial_grid$cell_state[i, j],
            free_virus1 = spatial_grid$free_virus1[i, j],
            ca_virus1 = spatial_grid$ca_virus1[i, j],
            free_virus2 = spatial_grid$free_virus2[i, j],
            ca_virus2 = spatial_grid$ca_virus2[i, j]
        )
        if(site$ca_virus2>0 && dice_cav2<=PARAMETERS$viral_ejection_rate*site$ca_virus2/PARAMETERS$max_rate){
          ca_viral_ejection(i,j,virus_type=2);
        }
      }
      
      
      ##Cell-free virus - virus 1
      ##--------------------------------------------
      #site<-spatial_grid[[key]];  #Re-fetch site after CA virus updates
      site <- list(
          cell_state = spatial_grid$cell_state[i, j],
          free_virus1 = spatial_grid$free_virus1[i, j],
          ca_virus1 = spatial_grid$ca_virus1[i, j],
          free_virus2 = spatial_grid$free_virus2[i, j],
          ca_virus2 = spatial_grid$ca_virus2[i, j]
      )
      if(site$free_virus1>0){
        
        #Event 0a: Decay
        viral_decay(i,j,dice_v1,virus_type=1);
        
        #Event 0b: Diffusion of remaining free virus at site
        if(dice_v1<=PARAMETERS$diff_rate_virus/PARAMETERS$max_rate){
          diffusion_virus(i,j,gaussian_mask_diffusion,virus_type=1);}			
        
      }
      
      ##Cell-free virus - virus 2
      ##--------------------------------------------
      #site<-spatial_grid[[key]];
      site <- list(
          cell_state = spatial_grid$cell_state[i, j],
          free_virus1 = spatial_grid$free_virus1[i, j],
          ca_virus1 = spatial_grid$ca_virus1[i, j],
          free_virus2 = spatial_grid$free_virus2[i, j],
          ca_virus2 = spatial_grid$ca_virus2[i, j]
      )
      if(site$free_virus2>0){
        
        #Event 0a: Decay
        viral_decay(i,j,dice_v2,virus_type=2);
        
        #Event 0b: Diffusion of remaining free virus at site
        if(dice_v2<=PARAMETERS$diff_rate_virus/PARAMETERS$max_rate){
          diffusion_virus(i,j,gaussian_mask_diffusion,virus_type=2);}			
        
      }
      
      
    }
    
    dynamic$num_updates<<-dynamic$num_updates+1;
    compute_totals();
    
    if(SETTINGS$visualize!='off'&&dynamic$num_updates==next_screen_grab){
      visualize();
      next_screen_grab<-dynamic$num_updates+SETTINGS$refr_freq;
    }
    output_results(results,F);
    if(dynamic$susceptible_cells==0 ||
       (dynamic$virions1==0&&dynamic$ca_virions1==0&&
        dynamic$virions2==0&&dynamic$ca_virions2==0&&
        dynamic$infected_nonprodcells==0&&dynamic$infected_prodcells==0)){
      #quit if there is no possibility of virus or targets exhausted
      break()
    }
  }
  
  end_simulation();
}

#PLOT DYNAMICS - creates summary plots from results.tsv
plot_dynamics<-function(){
  #Read results file
  results<-read.table(paste(SETTINGS$dirname,'/',SETTINGS$out_fname,sep=''),
                      header=TRUE,sep='\t',stringsAsFactors=FALSE);
  
  #Create plot file
  png(file=paste(SETTINGS$dirname,'/dynamics.png',sep=''),
      width=12,height=5,units='in',res=300);
  par(mfrow=c(1,3),mar=c(4,4,2,1));
  
  #Panel 1: Cell counts over time
  max_cells<-max(c(results$susceptible_cells,results$infected_nonprodcells,
                   results$infected_prodcells,results$dead_cells),na.rm=TRUE);
  plot(results$time,results$susceptible_cells,type='l',col='orange',lwd=2,
       xlab='Time (hours)',ylab='Cell count',main='Cell dynamics',
       ylim=c(0,max_cells));
  lines(results$time,results$infected_nonprodcells,col='turquoise',lwd=2);
  lines(results$time,results$infected_prodcells,col='blue',lwd=2);
  lines(results$time,results$dead_cells,col='black',lwd=2);
  legend('topright',
         legend=c('Susceptible','Infected (non-prod)','Infected (prod)','Dead'),
         col=c('orange','turquoise','blue','black'),
         lwd=2,bty='n',cex=0.8);
  
  #Panel 2: Virus counts over time (log scale) - separate by type
  max_virus<-max(c(results$virions1,results$ca_virions1,
                   results$virions2,results$ca_virions2),na.rm=TRUE);
  if(max_virus>0){
    plot(results$time,pmax(results$virions1+results$ca_virions1,1),type='l',col='darkgreen',lwd=2,
         xlab='Time (hours)',ylab='Virus count (log10)',main='Virus dynamics',
         log='y',ylim=c(1,max_virus*1.1));
    lines(results$time,pmax(results$virions2+results$ca_virions2,1),col='purple',lwd=2);
    legend('topright',
           legend=c('Virus 1 (total)','Virus 2 (total)'),
           col=c('darkgreen','purple'),
           lwd=2,bty='n',cex=0.8);
  }else{
    plot(results$time,rep(0,nrow(results)),type='l',
         xlab='Time (hours)',ylab='Virus count',main='Virus dynamics',
         ylim=c(0,1));
    text(mean(results$time),0.5,'No virus detected',cex=1.2);
  }
  
  #Panel 3: Coinfection over time
  if('coinfected_cells' %in% colnames(results)){
    plot(results$time,results$coinfected_cells,type='l',col='red',lwd=2,
         xlab='Time (hours)',ylab='Sites with both viruses',main='Coinfection');
  }
  
  dev.off();
  print(paste('Dynamics plot saved to',paste(SETTINGS$dirname,'/dynamics.png',sep='')));
}

#PLACE HOST cells
place_host<-function(){
  if(PARAMETERS$host_density==1){
    for(i in 1:SETTINGS$L){
      for(j in 1:SETTINGS$L){
        #key<-paste(i,j,sep=",");
        #spatial_grid[[key]]$cell_state<-1;
        spatial_grid$cell_state[i, j] <<- 1
      }
    }
  }else{
    for(i in 1:SETTINGS$L){
      for(j in 1:SETTINGS$L){
        if(runif(1)<PARAMETERS$host_density){
          #key<-paste(i,j,sep=",");
          #spatial_grid[[key]]$cell_state<-1;
          spatial_grid$cell_state[i, j] <<- 1
        }
      }
    }
  }
  gc();
}

#PLACE VIRUS
place_virus<-function(virus_type=1){
  place_settings<-if(virus_type==1) SETTINGS$place_virus1 else SETTINGS$place_virus2;
  
  if(place_settings$num==0){
    return(NULL); #do nothing and exit
  }else{
    if(place_settings$type=='single point'){ #place virus at the center of the grid
      center_i<-SETTINGS$L/2;
      center_j<-SETTINGS$L/2;
      #key<-paste(center_i,center_j,sep=",");
      if(virus_type==1){
        #spatial_grid[[key]]$free_virus1<-place_settings$num;
        spatial_grid$free_virus1[center_i, center_j] <<- place_settings$num;
      }else{
        #spatial_grid[[key]]$free_virus2<-place_settings$num;
        spatial_grid$free_virus2[center_i, center_j] <<- place_settings$num;
      }
      
    }else if(place_settings$type=='random'){ #virus spread randomly on grid
      for(v in 1:place_settings$num){
        i<-sample(1:SETTINGS$L,1);
        j<-sample(1:SETTINGS$L,1);
        #key<-paste(i,j,sep=",");
        if(virus_type==1){
          #spatial_grid[[key]]$free_virus1<-spatial_grid[[key]]$free_virus1+1;
          spatial_grid$free_virus1[i, j] <<- spatial_grid$free_virus1[i, j] + 1
        }else{
          #spatial_grid[[key]]$free_virus2<-spatial_grid[[key]]$free_virus2+1;
          spatial_grid$free_virus2[i, j] <<- spatial_grid$free_virus2[i, j] + 1
        }
      }
      
    }else if(place_settings$type=='infected producer'){ #place infected cell(s) instead of virus
      if(place_settings$num==1){ #at the center of the grid
        center_i<-SETTINGS$L/2;
        center_j<-SETTINGS$L/2;
        #key<-paste(center_i,center_j,sep=",");
        #spatial_grid[[key]]$cell_state<-3; 
        spatial_grid$cell_state[center_i, center_j] <<- 3
      }else{ #randomly on the grid
        for(v in 1:place_settings$num){
          i<-sample(1:SETTINGS$L,1);
          j<-sample(1:SETTINGS$L,1);
          #key<-paste(i,j,sep=",");
          #spatial_grid[[key]]$cell_state<-3;
          spatial_grid$cell_state[i, j] <<- 3
        }
      }
    }
  }
}

#Calculates plaque diameter in mm
calc_plaque_size<-function(){
  with(dynamic,
       return(sqrt(dead_cells)*SETTINGS$site_scale/1000));
}

#VISUALIZE
visualize<-function(){
  if(SETTINGS$visualize=='gif'){
    png(file=paste(SETTINGS$dirname,'/imgs/',
                   round(dynamic$time,digits=1),'.png',sep='')) ;
    par(mar=c(1,1,2,1),oma=c(1,0,0,0));
    layout(matrix(c(1:2),1,2,byrow = TRUE),widths=c(1,1),heights=c(1,1));
  }
  
  #Extract arrays from spatial_grid for visualization
  cell_state_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
  virus_type_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
  
  for(i in 1:SETTINGS$L){
    for(j in 1:SETTINGS$L){
      #key<-paste(i,j,sep=",");
      #cell_state_array[i,j]<-spatial_grid[[key]]$cell_state;
      cell_state_array[i, j] <- spatial_grid$cell_state[i, j]
      
      #Determine virus type at this site (including coinfection)
      #total_v1<-spatial_grid[[key]]$free_virus1+spatial_grid[[key]]$ca_virus1;
      total_v1 <- spatial_grid$free_virus1[i, j] + spatial_grid$ca_virus1[i, j]
      #total_v2<-spatial_grid[[key]]$free_virus2+spatial_grid[[key]]$ca_virus2;
      total_v2 <- spatial_grid$free_virus2[i, j] + spatial_grid$ca_virus2[i, j]
      
      if(total_v1==0 && total_v2==0){
        virus_type_array[i,j]<-0; #No virus
      }else if(total_v1>0 && total_v2>0){
        virus_type_array[i,j]<-3; #Coinfected
      }else if(total_v1>0){
        virus_type_array[i,j]<-1; #Virus 1 only
      }else{
        virus_type_array[i,j]<-2; #Virus 2 only
      }
    }
  }
  
  #Panel 1: Host cells 
  cols<-c('grey95','papayawhip','orange','turquoise','peachpuff4');  #0=Vacant, 1=uninf, 2=inf_np, 3=inf_p, 4=dead;
  types<-as.numeric(as.data.frame(table(cell_state_array),stringsAsFactors=FALSE)[,1]);
  image(cell_state_array,col=cols,breaks=c(-1:4),xaxt='n',yaxt='n',main='Cells')
  leg_text<-c('0: vacant','1: uninf','2: inf_np','3: inf_p','4: dead')
  legend('topright',leg_text[types+1],fill=cols[types+1],bty='n',ncol=1)
  # legend('bottomleft',paste(round(calc_plaque_size(settings),2),'mm'))
  
  #Panel 2: Virus by type
  virus_cols<-c('grey95','darkgreen','purple','red');
  if(sum(virus_type_array)==0){
    image(virus_type_array,col='grey95',xaxt='n',yaxt='n',main='Virus')
  }else{
    image(virus_type_array,col=virus_cols,breaks=c(-1:3),xaxt='n',yaxt='n',main='Virus by type');
    virus_types_present<-unique(as.vector(virus_type_array));
    virus_types_present<-virus_types_present[virus_types_present>0];
    if(length(virus_types_present)>0){
      leg_text<-c('Virus 1','Virus 2','Coinfected');
      legend('topright',leg_text[virus_types_present],
             fill=virus_cols[virus_types_present+1],bty='n',ncol=1);
    }
  }
  # legend('bottomleft',paste(dynamic$virions))
  
  #Timestamp
  mtext(paste('Time =',round(dynamic$time,digits=1),'h'),side=1,outer=T);
  
  if(SETTINGS$visualize=='gif'){dev.off()}
}

#COMPUTE GAUSSIAN kernel - only run once, at the beginning of simulation
compute_gaussian<-function(mask_range,sig=1.2,amp=1){
  temp<-(mask_range-1)/2;
  y<-matrix(rep(seq(-temp,temp),mask_range),mask_range); x<-t(y);
  gaussian_mask<-(amp/(2*pi*(sig^2)))*exp(-(x^2+y^2)/(2*sig^2));
  gaussian_mask<-gaussian_mask/sum(gaussian_mask); #this makes everything sum to 1
  return(gaussian_mask);
  rm(sig,y,x,temp); gc();
}

#DIFFUSION of free virus
diffusion_virus<-function(i,j,gaussian_mask,virus_type=1){
  #key<-paste(i,j,sep=",");
  #site<-spatial_grid[[key]];
  site <- list(
      cell_state = spatial_grid$cell_state[i, j],
      free_virus1 = spatial_grid$free_virus1[i, j],
      ca_virus1 = spatial_grid$ca_virus1[i, j],
      free_virus2 = spatial_grid$free_virus2[i, j],
      ca_virus2 = spatial_grid$ca_virus2[i, j]
  )
  
  #Determine which virus type to work with
  if(virus_type==1){
    virus_count<-site$free_virus1;
  }else{
    virus_count<-site$free_virus2;
  }
  
  if(virus_count==1){
    nbr<-pick_nbr(i,j);
    if(!is.null(nbr)){
      #nbr_key<-paste(nbr[1],nbr[2],sep=",");
      if(virus_type==1){
        #spatial_grid[[key]]$free_virus1<-0;
        spatial_grid$free_virus1[i, j] <<- 0
        #spatial_grid[[nbr_key]]$free_virus1<-spatial_grid[[nbr_key]]$free_virus1+1;
        spatial_grid$free_virus1[nbr[1], nbr[2]] <<- spatial_grid$free_virus1[nbr[1], nbr[2]] + 1
      }else{
        #spatial_grid[[key]]$free_virus2<-0;
        spatial_grid$free_virus2[i, j] <<- 0
        #spatial_grid[[nbr_key]]$free_virus2<-spatial_grid[[nbr_key]]$free_virus2+1;
        spatial_grid$free_virus2[nbr[1], nbr[2]] <<- spatial_grid$free_virus2[nbr[1], nbr[2]] + 1
      }
    }
  }else{
    burst<-gaussburst(virus_count,gaussian_mask);
    if(virus_type==1){
      #spatial_grid[[key]]$free_virus1<-0;
      spatial_grid$free_virus1[i, j] <<- 0
    }else{
      #spatial_grid[[key]]$free_virus2<-0;
      spatial_grid$free_virus2[i, j] <<- 0
    }
    distrib_progeny(i,j,burst,virus_type); 
  }
}

#GAUSSIAN BURST: generates a burst according to discrete approximation of a gaussian distribution
gaussburst<-function(burst_size,gaussian_mask){
  burst<-round(burst_size*gaussian_mask);
  if(sum(burst)>burst_size){ 
    exc<-sum(burst)-burst_size;
    tmp<-sample(seq(1,length(burst))[c(burst>0)],exc,replace=TRUE);
    counts <- tabulate(tmp, nbins = length(burst))
    burst <- burst - counts
    burst[burst<0]<-0;
  }else if(sum(burst)<burst_size){
    def<-burst_size-sum(burst);
    tmp<-sample(seq(1,length(burst)),def,replace=TRUE);
    counts <- tabulate(tmp, nbins = length(burst))
    burst <- burst + counts
  }
  return(burst);
}

#DISTRIBUTE PROGENY in free virus array for burst or diffusion
distrib_progeny <- function(i, j, burst, virus_type=1) {
    
    # determine the splatter zone!
    dist <- (PARAMETERS$diff_range - 1) / 2
    
    # bounds for the cells we're updating
    imin <- max(i - dist, 1)
    imax <- min(i + dist, SETTINGS$L)
    jmin <- max(j - dist, 1)
    jmax <- min(j + dist, SETTINGS$L)
    
    # iteratize over the cells to update
    for (a in seq(from = imin, to = imax, by = 1)) {
        for (b in seq(from = jmin, to = jmax, by = 1)) {
            
            # indices for burst matrix
            burst_i <- a - i + (dist + 1)
            burst_j <- b - j + (dist + 1)
            
            # increment the value of free_virus1 or 2 by the corresponding value in burst
            if(virus_type == 1) {
                spatial_grid$free_virus1[a, b] <<- spatial_grid$free_virus1[a, b] + burst[burst_i, burst_j]
            } else {
                spatial_grid$free_virus2[a, b] <<- spatial_grid$free_virus2[a, b] + burst[burst_i, burst_j]
            }
        }
    }
}

#PICK NEIGHBOUR from immediate moore neighbourhood
pick_nbr<-function(i,j){
  delta<-sample(c(-1,0,1),2,replace=TRUE);
  nbr<-c(i,j)+delta;
  if(any(nbr<1|nbr>SETTINGS$L)||sum(delta)==0){
    return(NULL); #outside grid
  }else{
    return(nbr);
  }
}

#DECAY of free virus
viral_decay<-function(i,j,dice,virus_type=1){
  #key<-paste(i,j,sep=",");
  site <- list(
      cell_state = spatial_grid$cell_state[i, j],
      free_virus1 = spatial_grid$free_virus1[i, j],
      ca_virus1 = spatial_grid$ca_virus1[i, j],
      free_virus2 = spatial_grid$free_virus2[i, j],
      ca_virus2 = spatial_grid$ca_virus2[i, j]
  )
  
  if(virus_type==1){
    if(site$free_virus1==1){
      if(dice<PARAMETERS$freevirus_decay/PARAMETERS$max_rate){
        #spatial_grid[[key]]$free_virus1<-0;
        spatial_grid$free_virus1[i, j] <<- 0
      }
    }else{
      #spatial_grid[[key]]$free_virus1<-round(site$free_virus1*exp(-PARAMETERS$freevirus_decay/PARAMETERS$max_rate));
      spatial_grid$free_virus1[i, j] <<- round(site$free_virus1*exp(-PARAMETERS$freevirus_decay/PARAMETERS$max_rate));
    }
  }else{
    if(site$free_virus2==1){
      if(dice<PARAMETERS$freevirus_decay/PARAMETERS$max_rate){
        #spatial_grid[[key]]$free_virus2<-0;
        spatial_grid$free_virus2[i, j] <<- 0
      }
    }else{
      #spatial_grid[[key]]$free_virus2<-round(site$free_virus2*exp(-PARAMETERS$freevirus_decay/PARAMETERS$max_rate));
      spatial_grid$free_virus2[i, j] <<- round(site$free_virus2*exp(-PARAMETERS$freevirus_decay/PARAMETERS$max_rate));
    }
  }
}

#DECAY of cell-associated virus
ca_viral_decay<-function(i,j,dice,virus_type=1){
  #key<-paste(i,j,sep=",");
  #site<-spatial_grid[[key]];
  site <- list(
      cell_state = spatial_grid$cell_state[i, j],
      free_virus1 = spatial_grid$free_virus1[i, j],
      ca_virus1 = spatial_grid$ca_virus1[i, j],
      free_virus2 = spatial_grid$free_virus2[i, j],
      ca_virus2 = spatial_grid$ca_virus2[i, j]
  )
  
  if(virus_type==1){
    if(site$ca_virus1==1){
      if(dice<PARAMETERS$cavirus_decay/PARAMETERS$max_rate){
        #spatial_grid[[key]]$ca_virus1<-0;
        spatial_grid$ca_virus1[i, j] <<- 0
      }
    }else{
      #spatial_grid[[key]]$ca_virus1<-round(site$ca_virus1*exp(-PARAMETERS$cavirus_decay/PARAMETERS$max_rate));
      spatial_grid$ca_virus1[i, j] <<- round(site$ca_virus1*exp(-PARAMETERS$cavirus_decay/PARAMETERS$max_rate));
    }
  }else{
    if(site$ca_virus2==1){
      if(dice<PARAMETERS$cavirus_decay/PARAMETERS$max_rate){
        #spatial_grid[[key]]$ca_virus2<-0;
        spatial_grid$ca_virus2[i, j] <<- 0
      }
    }else{
      #spatial_grid[[key]]$ca_virus2<-round(site$ca_virus2*exp(-PARAMETERS$cavirus_decay/PARAMETERS$max_rate));
      spatial_grid$ca_virus2[i, j] <<- round(site$ca_virus2*exp(-PARAMETERS$cavirus_decay/PARAMETERS$max_rate));
    }
  }
}

#EJECTION of cell-associated virus to free virus
ca_viral_ejection<-function(i,j,virus_type=1){
  # key<-paste(i,j,sep=",");
  # site<-spatial_grid[[key]];
  site <- list(
      cell_state = spatial_grid$cell_state[i, j],
      free_virus1 = spatial_grid$free_virus1[i, j],
      ca_virus1 = spatial_grid$ca_virus1[i, j],
      free_virus2 = spatial_grid$free_virus2[i, j],
      ca_virus2 = spatial_grid$ca_virus2[i, j]
  )
  
  if(virus_type==1 && site$ca_virus1>0){
    num_to_eject<-round(site$ca_virus1*PARAMETERS$ca_to_free_proportion);
    if(num_to_eject>0){
      #spatial_grid[[key]]$ca_virus1<-site$ca_virus1-num_to_eject;
      spatial_grid$ca_virus1[i, j] <<- site$ca_virus1 - num_to_eject
      #spatial_grid[[key]]$free_virus1<-site$free_virus1+num_to_eject;
      spatial_grid$free_virus1[i, j] <<- site$free_virus1 + num_to_eject
    }
  }else if(virus_type==2 && site$ca_virus2>0){
    num_to_eject<-round(site$ca_virus2*PARAMETERS$ca_to_free_proportion);
    if(num_to_eject>0){
      #spatial_grid[[key]]$ca_virus2<-site$ca_virus2-num_to_eject;
      spatial_grid$ca_virus2[i, j] <<- site$ca_virus2 - num_to_eject
      #spatial_grid[[key]]$free_virus2<-site$free_virus2+num_to_eject;
      spatial_grid$free_virus2[i, j] <<- site$free_virus2 + num_to_eject
    }
  }
}

# OUTPUT 
output_results<-function(results,header){
  output_file<-paste(SETTINGS$dirname,'/',SETTINGS$out_fname,sep='');
  
  if (!nchar(output_file) || output_file == "") {
    stop('No input file');
  }
  if (header) {
    cat(file=output_file,paste(names(dynamic),collapse='\t'),'\n',append=F);
  } else {
    cat(file=output_file,paste(dynamic,collapse='\t'),'\n',append=T);
  }
}

# Compute totals of cells and virus
compute_totals<-function(){
  
  dynamic$time<<-dynamic$num_updates/PARAMETERS$max_rate;
  
  dynamic$virions1    <<- sum(spatial_grid$free_virus1)
  dynamic$ca_virions1 <<- sum(spatial_grid$ca_virus1)
  dynamic$virions2    <<- sum(spatial_grid$free_virus2)
  dynamic$ca_virions2 <<- sum(spatial_grid$ca_virus2)
  
  dynamic$viral_cells1    <<- sum(spatial_grid$free_virus1 > 0)
  dynamic$ca_viral_cells1 <<- sum(spatial_grid$ca_virus1   > 0)
  dynamic$viral_cells2    <<- sum(spatial_grid$free_virus2 > 0)
  dynamic$ca_viral_cells2 <<- sum(spatial_grid$ca_virus2   > 0)
  
  # coinfections (sites with both virus types present (free or CA))
  dynamic$coinfected_cells <<- sum(
      (spatial_grid$free_virus1 > 0 | spatial_grid$ca_virus1 > 0) &
      (spatial_grid$free_virus2 > 0 | spatial_grid$ca_virus2 > 0)
  )
  
  dynamic$susceptible_cells     <<- sum(spatial_grid$cell_state == 1)
  dynamic$infected_nonprodcells <<- sum(spatial_grid$cell_state == 2)
  dynamic$infected_prodcells    <<- sum(spatial_grid$cell_state == 3)
  dynamic$dead_cells            <<- sum(spatial_grid$cell_state == 4)
  
  total_virions<-dynamic$virions1+dynamic$virions2;
  if (total_virions > dynamic$max_vl) {
    dynamic$max_vl <<-  total_virions;
    dynamic$max_vl_time <<-  dynamic$time;
  }
}

#END of simulation - save grids
end_simulation<-function(){
  if(SETTINGS$save_grids){
    #Extract arrays from spatial_grid structure for saving
    cell_state_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
    free_virus1_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
    ca_virus1_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
    free_virus2_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
    ca_virus2_array<-matrix(0,nrow=SETTINGS$L,ncol=SETTINGS$L);
    
    for(i in 1:SETTINGS$L){
      for(j in 1:SETTINGS$L){
        # key<-paste(i,j,sep=",");
        # cell_state_array[i,j]<-spatial_grid[[key]]$cell_state;
        cell_state_array[i, j] <- spatial_grid$cell_state[i, j]
        # free_virus1_array[i,j]<-spatial_grid[[key]]$free_virus1;
        free_virus1_array[i, j] <- spatial_grid$free_virus1[i, j]
        # ca_virus1_array[i,j]<-spatial_grid[[key]]$ca_virus1;
        ca_virus1_array[i, j] <- spatial_grid$ca_virus1[i, j]
        # free_virus2_array[i,j]<-spatial_grid[[key]]$free_virus2;
        free_virus2_array[i, j] <- spatial_grid$free_virus2[i, j]
        # ca_virus2_array[i,j]<-spatial_grid[[key]]$ca_virus2;
        ca_virus2_array[i, j] <- spatial_grid$ca_virus2[i, j]
      }
    }
    save(cell_state_array,file=paste(SETTINGS$dirname,'/cell_state',sep=''));
    save(free_virus1_array,file=paste(SETTINGS$dirname,'/free_virus1',sep=''));
    save(ca_virus1_array,file=paste(SETTINGS$dirname,'/ca_virus1',sep=''));
    save(free_virus2_array,file=paste(SETTINGS$dirname,'/free_virus2',sep=''));
    save(ca_virus2_array,file=paste(SETTINGS$dirname,'/ca_virus2',sep=''));
  }
  if(SETTINGS$visualize=='pdf'){
    dev.off();
  }
  
  #Plot dynamics
  plot_dynamics();
  
  print(paste("Highest log VL",ifelse(dynamic$max_vl>0,log10(dynamic$max_vl),0),"at t=",dynamic$max_vl_time));
  print(paste("End of simulation at t=",dynamic$time));
}

# run ---------------------------------------------------------------------
set.seed(123)

library(profvis)
profvis({
    Spatial_main()
})

