# Spatial model of HSV coinfection in skin (2D) - List-based version
# Pavitra Roychoudhury
# Fred Hutchinson Cancer Research Center
# Tracks two virus types and coinfection without recombination

# INITIALIZATION: lines to run when this file is sourced
rm(list=ls()); 
library(RColorBrewer);

# 0.  MAIN FUNCTION-only used if there is no wrapper file
Spatial_main<-function(visualize='off',param_file=NULL,settings_file=NULL){
  options(error=recover);
  settings<-get_default_settings();
  parameters<-get_default_parameters();
  
  #Load settings from file if provided
  if(!is.null(settings_file)){
    settings<-load_settings_from_file(settings_file,settings);
  }
  
  #Load parameters from file if provided
  if(!is.null(param_file)){
    parameters<-load_parameters_from_file(param_file,parameters);
  }
  
  #modify settings and parameters here or use a wrapper
  settings$visualize<-visualize;
  
  #initialize and run sim
  settings<-init_grid(settings,parameters);
  simulation(settings,parameters);
}

## FUNCTIONS: all the functions that are needed for the simulation are below this line

# 1a.  Get settings for the simulation
get_default_settings<-function(){
  settings<-list(L=200);                          #Specify grid dimensions
  settings$dirname<-'Results/';   								#Directory to store results
  settings$out_fname<-'results.tsv';							#File to store cell and virus counts
  settings$sim_length<-1;											#Length of simulation in hours
  settings$save_grids<-TRUE;          #Save grids at the end of the sim? 
  
  #Visualization
  settings$visualize<-'off';                         #Options: 'screen','pdf','off','gif'
  settings$refr_freq<-25;												#Refresh frequency (hrs)
  settings$site_scale<-50; #how many microns is one site
  
  #Place virus
  settings$place_virus1<-list(type='infected producer',num=10);  #Initial virus 1
  settings$place_virus2<-list(type='infected producer',num=10);  #Initial virus 2
  
  return(settings);
}

# 1b2. Load settings from CSV file
# CSV file should have two columns: setting_name, value
# Example:
# setting_name,value
# L,150
# sim_length,10
# save_grids,TRUE
# refr_freq,50
load_settings_from_file<-function(settings_file,settings){
  if(!file.exists(settings_file)){
    stop(paste('Settings file not found:',settings_file));
  }
  
  #Read CSV file
  settings_data<-read.csv(settings_file,stringsAsFactors=FALSE);
  
  #Check that file has correct columns
  if(!all(c('setting_name','value') %in% colnames(settings_data))){
    stop('Settings file must have columns: setting_name, value');
  }
  
  #Update settings
  for(i in 1:nrow(settings_data)){
    setting_name<-settings_data$setting_name[i];
    setting_value<-settings_data$value[i];
    
    #Check if setting exists
    if(setting_name %in% names(settings)){
      #Convert to appropriate type
      if(is.numeric(settings[[setting_name]])){
        setting_value<-as.numeric(setting_value);
      }else if(is.logical(settings[[setting_name]])){
        setting_value<-as.logical(setting_value);
      }
      settings[[setting_name]]<-setting_value;
      print(paste('Updated setting:',setting_name,'=',setting_value));
    }else{
      warning(paste('Setting',setting_name,'not found in default settings. Skipping.'));
    }
  }
  
  return(settings);
}

# 1c. Get default parameters 
get_default_parameters<-function(){
  parameters<-list(); 
  
  #Rates (per site or cell per day)
  parameters$max_rate<-25;						  #Sets time resolution, should be greater than or equal to sum of all rates for a given state
  
  #Uninfected cells
  parameters$uninf_cell_death<-0;	#Rate of uninfected cell death 
  parameters$infectivity_free<-0.1;				#Infectivity (per virus per cell per day), will get multiplied by number of free viruses at the site
  parameters$cell_div<-0.077; 					#Rate of generation of new susceptible cells (repop of empty sites)
  
  #Infected cells
  parameters$inf_cell_death<-1.25;			#Rate of infected cell death due to infection
  parameters$viral_lag<-8;						  #Rate at which inf non-producer becomes viral producing cell (per cell per day)
  parameters$vir_prod_rate_mean<-100000;	#Rate at which an infected cell produces free virus (per cell per day, mean and sd)
  parameters$vir_prod_rate_sd<-10000;
  
  #Virus
  parameters$diff_rate_virus<-12; 						  #Rate of diffusion of free virus (sites moved per day, depends on diffusivity, size of virion, viscosity of medium, etc)
  parameters$freevirus_decay<-8.8;		#Rate of decay of free virus (per day)
  parameters$cavirus_decay<-8.8;		#Rate of decay of cell-associated virus (per day)
  parameters$viral_ejection_rate<-10;  #Rate of conversion from cell-associated to free virus (per day)
  parameters$ca_to_free_proportion<-0.9;  #Proportion of CA virus converted to free virus per ejection event
  
  #Probabilities, fractions, etc
  parameters$host_density<-1;           #Fraction of sites occupied by susceptible cells at start	
  parameters$diff_range<-3;						  #n x n site, default n = 3, i.e. immediate Moore neighbourhood. 
  
  return(parameters);
}

# 1d. Load parameters from CSV file
# CSV file should have two columns: parameter_name, value
# Example:
# parameter_name,value
# max_rate,30
# infectivity_free,0.15
load_parameters_from_file<-function(param_file,parameters){
  if(!file.exists(param_file)){
    stop(paste('Parameter file not found:',param_file));
  }
  
  #Read CSV file
  param_data<-read.csv(param_file,stringsAsFactors=FALSE);
  
  #Check that file has correct columns
  if(!all(c('parameter_name','value') %in% colnames(param_data))){
    stop('Parameter file must have columns: parameter_name, value');
  }
  
  #Update parameters
  for(i in 1:nrow(param_data)){
    param_name<-param_data$parameter_name[i];
    param_value<-param_data$value[i];
    
    #Check if parameter exists
    if(param_name %in% names(parameters)){
      #Convert to numeric if needed
      if(is.numeric(parameters[[param_name]])){
        param_value<-as.numeric(param_value);
      }
      parameters[[param_name]]<-param_value;
      print(paste('Updated parameter:',param_name,'=',param_value));
    }else{
      warning(paste('Parameter',param_name,'not found in default parameters. Skipping.'));
    }
  }
  
  return(parameters);
}

# 2. Initialize grid and other arrays
init_grid<-function(settings,parameters){
  print('Starting simulation.. hold on to your hat!');
  
  #Compute simulation length in number of updates
  settings$tot_updates<-settings$sim_length*parameters$max_rate;
  
  #Create spatial grid as an environment for fast hashed access
  #Each grid location is accessed via key "i,j" and contains a list with cell_state and virus counts for each type
  spatial_grid<<-new.env(hash=TRUE);
  for(i in 1:settings$L){
    for(j in 1:settings$L){
      key<-paste(i,j,sep=",");
      spatial_grid[[key]]<<-list(cell_state=as.integer(0), 
                                 free_virus1=as.integer(0), ca_virus1=as.integer(0),
                                 free_virus2=as.integer(0), ca_virus2=as.integer(0));
    }
  }
  
  #Directories
  print(paste('Results directory =',settings$dirname));
  dir.create(settings$dirname,showWarnings=FALSE);
  save(parameters,file=paste(settings$dirname,'/simulationparameters',sep=''));
  save(settings,file=paste(settings$dirname,'/simulationsettings',sep=''));
  if(settings$visualize=='gif')dir.create(paste(settings$dirname,'/imgs',sep=''),showWarnings=FALSE)
  
  #Dynamic variables - this stores time elapsed and population counts
  dynamic<<-list(num_updates=0,time=0,
                 susceptible_cells=0,infected_nonprodcells=0,infected_prodcells=0,dead_cells=0,
                 virions1=0,ca_virions1=0,virions2=0,ca_virions2=0,
                 max_vl=0,max_vl_time=0,
                 viral_cells1=0,ca_viral_cells1=0,viral_cells2=0,ca_viral_cells2=0,
                 coinfected_cells=0);
  output_results(settings,dynamic,T); #write header
  
  #Visualization setup
  if(settings$visualize!='off'){
    if(settings$visualize=='pdf'){
      pdf(file=paste(settings$dirname,'/Results.pdf',sep=''),
          width=8,height=4); 
    }else if(settings$visualize=='screen'){
      quartz(width=8,height=4);   
    }else if(settings$visualize=='gif'){
      # library(animation)
      # png(file=paste(settings$dirname,'/imgs/Results.png',sep=''))
    }else{
      stop('Not a valid setting for "visualize"')
    }   	
    par(mar=c(1,1,2,1),oma=c(1,0,0,0));
    layout(matrix(c(1:2),1,2,byrow = TRUE),widths=c(1,1),heights=c(1,1));
  }
  gc(); return(settings)
}

# 3. SIMULATION
simulation<-function(settings,parameters){
  
  gaussian_mask_diffusion<-compute_gaussian(parameters$diff_range); 
  
  place_host(settings,parameters); #Place host cells
  if(settings$place_virus1$type!='none') place_virus(settings,virus_type=1); #Place virus 1
  if(settings$place_virus2$type!='none') place_virus(settings,virus_type=2); #Place virus 2
  
  dynamic$num_updates<<-0;
  compute_totals(settings,parameters);
  output_results(settings,dynamic,F);
  if(settings$visualize!='off'){
    visualize(settings); 
  }
  next_screen_grab<-dynamic$num_updates+settings$refr_freq;
  
  while(dynamic$num_updates<settings$tot_updates){
    
    #Generate random order of sites to update
    update_order<-sample(seq(1,settings$L^2),settings$L^2,replace=TRUE)
    row_order<-((update_order-1)%%settings$L)+1;
    col_order<-floor((update_order-1)/settings$L)+1;
    dice_cells_all<-runif(settings$L^2,min=0,max=1);
    dice_v1_all<-runif(settings$L^2,min=0,max=1);
    dice_v2_all<-runif(settings$L^2,min=0,max=1);
    dice_cav1_all<-runif(settings$L^2,min=0,max=1);
    dice_cav2_all<-runif(settings$L^2,min=0,max=1);
    
    #Update sites in random order (same as picking a random site to update one at a time)
    for(rep in 1:settings$L^2){
      i<-row_order[rep];
      j<-col_order[rep];
      dice_cells<-dice_cells_all[rep];
      dice_v1<-dice_v1_all[rep];
      dice_v2<-dice_v2_all[rep];
      dice_cav1<-dice_cav1_all[rep];
      dice_cav2<-dice_cav2_all[rep];
      
      key<-paste(i,j,sep=",");
      site<-spatial_grid[[key]];
      
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
        if(dice_cells<=parameters$cell_div/parameters$max_rate){
          spatial_grid[[key]]$cell_state<-1;
        }
        
        
        #State 1: uninfected cell 
      }else if(site$cell_state==1){
        
        #Event 1.1: infection due to free virus at the site (any virus type can infect)
        total_free_virus<-site$free_virus1+site$free_virus2;
        if(total_free_virus>0 && dice_cells<=parameters$infectivity_free*total_free_virus/parameters$max_rate){
          spatial_grid[[key]]$cell_state<-2; #becomes an infected cell
          #Remove one virus of random type
          virus_probs<-c(site$free_virus1,site$free_virus2)/total_free_virus;
          infecting_type<-sample(1:2,1,prob=virus_probs);
          if(infecting_type==1){
            spatial_grid[[key]]$free_virus1<-site$free_virus1-1;
          }else{
            spatial_grid[[key]]$free_virus2<-site$free_virus2-1;
          }
          
          #Event 1.2: apoptosis
        }else if(dice_cells<=parameters$uninf_cell_death/parameters$max_rate){
          spatial_grid[[key]]$cell_state<-4;
        } 
        
        #State 2: infected cell (non-virus producing)
      }else if(site$cell_state==2){					
        
        #Event 2.1: become a virus-producing cell
        if(dice_cells<=parameters$viral_lag/parameters$max_rate){
          spatial_grid[[key]]$cell_state<-3; #becomes virus-producing infected cell
          
          #Event 2.2: dies 
        }else if(dice_cells<=(parameters$viral_lag+parameters$uninf_cell_death)/parameters$max_rate){
          spatial_grid[[key]]$cell_state<-4; #death of infected cell
        }
        
        #State 3: infected cell (producing virus)
      }else if(site$cell_state==3){
        
        #Event 3.1: die
        if(dice_cells<=parameters$inf_cell_death/parameters$max_rate){
          spatial_grid[[key]]$cell_state<-4; #death of infected cell
          
          #Event 3.2: produce cell-associated virus (virus 1)
        }else{
          new_virus1<-max(0,round(rnorm(1,mean=parameters$vir_prod_rate_mean/parameters$max_rate,
                                        sd=parameters$vir_prod_rate_sd/parameters$max_rate)));
          spatial_grid[[key]]$ca_virus1<-site$ca_virus1+new_virus1;
          
          #Event 3.3: produce cell-associated virus (virus 2)
          new_virus2<-max(0,round(rnorm(1,mean=parameters$vir_prod_rate_mean/parameters$max_rate,
                                        sd=parameters$vir_prod_rate_sd/parameters$max_rate)));
          spatial_grid[[key]]$ca_virus2<-site$ca_virus2+new_virus2;
        }
      }
      
      
      ##Cell-associated virus - virus 1
      ##--------------------------------------------
      site<-spatial_grid[[key]];  #Re-fetch site after cell state updates
      if(site$ca_virus1>0){
        
        #Event 0a: Decay
        ca_viral_decay(i,j,parameters,dice_cav1,virus_type=1);
        
        #Event 0b: Ejection to free virus (probability increases with number of CA viruses)
        site<-spatial_grid[[key]];  #Re-fetch after decay
        if(site$ca_virus1>0 && dice_cav1<=parameters$viral_ejection_rate*site$ca_virus1/parameters$max_rate){
          ca_viral_ejection(i,j,parameters,virus_type=1);
        }
      }
      
      ##Cell-associated virus - virus 2
      ##--------------------------------------------
      site<-spatial_grid[[key]];
      if(site$ca_virus2>0){
        
        #Event 0a: Decay
        ca_viral_decay(i,j,parameters,dice_cav2,virus_type=2);
        
        #Event 0b: Ejection to free virus
        site<-spatial_grid[[key]];
        if(site$ca_virus2>0 && dice_cav2<=parameters$viral_ejection_rate*site$ca_virus2/parameters$max_rate){
          ca_viral_ejection(i,j,parameters,virus_type=2);
        }
      }
      
      
      ##Cell-free virus - virus 1
      ##--------------------------------------------
      site<-spatial_grid[[key]];  #Re-fetch site after CA virus updates
      if(site$free_virus1>0){
        
        #Event 0a: Decay
        viral_decay(i,j,parameters,dice_v1,virus_type=1);
        
        #Event 0b: Diffusion of remaining free virus at site
        if(dice_v1<=parameters$diff_rate_virus/parameters$max_rate){
          diffusion_virus(settings,i,j,gaussian_mask_diffusion,parameters,virus_type=1);}			
        
      }
      
      ##Cell-free virus - virus 2
      ##--------------------------------------------
      site<-spatial_grid[[key]];
      if(site$free_virus2>0){
        
        #Event 0a: Decay
        viral_decay(i,j,parameters,dice_v2,virus_type=2);
        
        #Event 0b: Diffusion of remaining free virus at site
        if(dice_v2<=parameters$diff_rate_virus/parameters$max_rate){
          diffusion_virus(settings,i,j,gaussian_mask_diffusion,parameters,virus_type=2);}			
        
      }
      
      
    }
    
    dynamic$num_updates<<-dynamic$num_updates+1;
    compute_totals(settings,parameters);
    
    if(settings$visualize!='off'&&dynamic$num_updates==next_screen_grab){
      visualize(settings);
      next_screen_grab<-dynamic$num_updates+settings$refr_freq;
    }
    output_results(settings,results,F);
    if(dynamic$susceptible_cells==0 ||
       (dynamic$virions1==0&&dynamic$ca_virions1==0&&
        dynamic$virions2==0&&dynamic$ca_virions2==0&&
        dynamic$infected_nonprodcells==0&&dynamic$infected_prodcells==0)){
      #quit if there is no possibility of virus or targets exhausted
      break()
    }
  }
  
  end_simulation(settings,parameters);
}

#PLOT DYNAMICS - creates summary plots from results.tsv
plot_dynamics<-function(settings){
  #Read results file
  results<-read.table(paste(settings$dirname,'/',settings$out_fname,sep=''),
                      header=TRUE,sep='\t',stringsAsFactors=FALSE);
  
  #Create plot file
  png(file=paste(settings$dirname,'/dynamics.png',sep=''),
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
  print(paste('Dynamics plot saved to',paste(settings$dirname,'/dynamics.png',sep='')));
}

#PLACE HOST cells
place_host<-function(settings,parameters){
  if(parameters$host_density==1){
    for(i in 1:settings$L){
      for(j in 1:settings$L){
        key<-paste(i,j,sep=",");
        spatial_grid[[key]]$cell_state<-1;
      }
    }
  }else{
    for(i in 1:settings$L){
      for(j in 1:settings$L){
        if(runif(1)<parameters$host_density){
          key<-paste(i,j,sep=",");
          spatial_grid[[key]]$cell_state<-1;
        }
      }
    }
  }
  gc();
}

#PLACE VIRUS
place_virus<-function(settings,virus_type=1){
  place_settings<-if(virus_type==1) settings$place_virus1 else settings$place_virus2;
  
  if(place_settings$num==0){
    return(NULL); #do nothing and exit
  }else{
    if(place_settings$type=='single point'){ #place virus at the center of the grid
      center_i<-settings$L/2;
      center_j<-settings$L/2;
      key<-paste(center_i,center_j,sep=",");
      if(virus_type==1){
        spatial_grid[[key]]$free_virus1<-place_settings$num;
      }else{
        spatial_grid[[key]]$free_virus2<-place_settings$num;
      }
      
    }else if(place_settings$type=='random'){ #virus spread randomly on grid
      for(v in 1:place_settings$num){
        i<-sample(1:settings$L,1);
        j<-sample(1:settings$L,1);
        key<-paste(i,j,sep=",");
        if(virus_type==1){
          spatial_grid[[key]]$free_virus1<-spatial_grid[[key]]$free_virus1+1;
        }else{
          spatial_grid[[key]]$free_virus2<-spatial_grid[[key]]$free_virus2+1;
        }
      }
      
    }else if(place_settings$type=='infected producer'){ #place infected cell(s) instead of virus
      if(place_settings$num==1){ #at the center of the grid
        center_i<-settings$L/2;
        center_j<-settings$L/2;
        key<-paste(center_i,center_j,sep=",");
        spatial_grid[[key]]$cell_state<-3; 
      }else{ #randomly on the grid
        for(v in 1:place_settings$num){
          i<-sample(1:settings$L,1);
          j<-sample(1:settings$L,1);
          key<-paste(i,j,sep=",");
          spatial_grid[[key]]$cell_state<-3;
        }
      }
    }
  }
}

#Calculates plaque diameter in mm
calc_plaque_size<-function(settings){
  with(dynamic,
       return(sqrt(dead_cells)*settings$site_scale/1000));
}

#VISUALIZE
visualize<-function(settings){
  if(settings$visualize=='gif'){
    png(file=paste(settings$dirname,'/imgs/',
                   round(dynamic$time,digits=1),'.png',sep='')) ;
    par(mar=c(1,1,2,1),oma=c(1,0,0,0));
    layout(matrix(c(1:2),1,2,byrow = TRUE),widths=c(1,1),heights=c(1,1));
  }
  
  #Extract arrays from spatial_grid for visualization
  cell_state_array<-matrix(0,nrow=settings$L,ncol=settings$L);
  virus_type_array<-matrix(0,nrow=settings$L,ncol=settings$L);
  
  for(i in 1:settings$L){
    for(j in 1:settings$L){
      key<-paste(i,j,sep=",");
      cell_state_array[i,j]<-spatial_grid[[key]]$cell_state;
      
      #Determine virus type at this site (including coinfection)
      total_v1<-spatial_grid[[key]]$free_virus1+spatial_grid[[key]]$ca_virus1;
      total_v2<-spatial_grid[[key]]$free_virus2+spatial_grid[[key]]$ca_virus2;
      
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
  
  if(settings$visualize=='gif'){dev.off()}
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
diffusion_virus<-function(settings,i,j,gaussian_mask,parameters,virus_type=1){
  key<-paste(i,j,sep=",");
  site<-spatial_grid[[key]];
  
  #Determine which virus type to work with
  if(virus_type==1){
    virus_count<-site$free_virus1;
  }else{
    virus_count<-site$free_virus2;
  }
  
  if(virus_count==1){
    nbr<-pick_nbr(settings,i,j);
    if(!is.null(nbr)){
      nbr_key<-paste(nbr[1],nbr[2],sep=",");
      if(virus_type==1){
        spatial_grid[[key]]$free_virus1<-0;
        spatial_grid[[nbr_key]]$free_virus1<-spatial_grid[[nbr_key]]$free_virus1+1;
      }else{
        spatial_grid[[key]]$free_virus2<-0;
        spatial_grid[[nbr_key]]$free_virus2<-spatial_grid[[nbr_key]]$free_virus2+1;
      }
    }
  }else{
    burst<-gaussburst(virus_count,gaussian_mask);
    if(virus_type==1){
      spatial_grid[[key]]$free_virus1<-0;
    }else{
      spatial_grid[[key]]$free_virus2<-0;
    }
    distrib_progeny(settings,i,j,burst,parameters$diff_range,parameters,virus_type); 
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
distrib_progeny<-function(settings,i,j,burst,dist_range,parameters,virus_type=1){
  temp2<-(dist_range-1)/2;
  imin<-max(1,i-temp2);imax<-min(i+temp2,settings$L); #deal with edges
  jmin<-max(1,j-temp2);jmax<-min(j+temp2,settings$L);
  imin1<-temp2-(i-imin)+1;imax1<-temp2+(imax-i)+1;
  jmin1<-temp2-(j-jmin)+1;jmax1<-temp2+(jmax-j)+1;
  
  for(ii in 1:length(imin:imax)){
    for(jj in 1:length(jmin:jmax)){
      actual_i<-imin+ii-1;
      actual_j<-jmin+jj-1;
      burst_i<-imin1+ii-1;
      burst_j<-jmin1+jj-1;
      key<-paste(actual_i,actual_j,sep=",");
      if(virus_type==1){
        spatial_grid[[key]]$free_virus1<-spatial_grid[[key]]$free_virus1+burst[burst_i,burst_j];
      }else{
        spatial_grid[[key]]$free_virus2<-spatial_grid[[key]]$free_virus2+burst[burst_i,burst_j];
      }
    }
  }
}

#PICK NEIGHBOUR from immediate moore neighbourhood
pick_nbr<-function(settings,i,j){
  delta<-sample(c(-1,0,1),2,replace=TRUE);
  nbr<-c(i,j)+delta;
  if(any(nbr<1|nbr>settings$L)||sum(delta)==0){
    return(NULL); #outside grid
  }else{
    return(nbr);
  }
}

#DECAY of free virus
viral_decay<-function(i,j,parameters,dice,virus_type=1){
  key<-paste(i,j,sep=",");
  site<-spatial_grid[[key]];
  
  if(virus_type==1){
    if(site$free_virus1==1){
      if(dice<parameters$freevirus_decay/parameters$max_rate){
        spatial_grid[[key]]$free_virus1<-0;
      }
    }else{
      spatial_grid[[key]]$free_virus1<-round(site$free_virus1*exp(-parameters$freevirus_decay/parameters$max_rate));
    }
  }else{
    if(site$free_virus2==1){
      if(dice<parameters$freevirus_decay/parameters$max_rate){
        spatial_grid[[key]]$free_virus2<-0;
      }
    }else{
      spatial_grid[[key]]$free_virus2<-round(site$free_virus2*exp(-parameters$freevirus_decay/parameters$max_rate));
    }
  }
}

#DECAY of cell-associated virus
ca_viral_decay<-function(i,j,parameters,dice,virus_type=1){
  key<-paste(i,j,sep=",");
  site<-spatial_grid[[key]];
  
  if(virus_type==1){
    if(site$ca_virus1==1){
      if(dice<parameters$cavirus_decay/parameters$max_rate){
        spatial_grid[[key]]$ca_virus1<-0;
      }
    }else{
      spatial_grid[[key]]$ca_virus1<-round(site$ca_virus1*exp(-parameters$cavirus_decay/parameters$max_rate));
    }
  }else{
    if(site$ca_virus2==1){
      if(dice<parameters$cavirus_decay/parameters$max_rate){
        spatial_grid[[key]]$ca_virus2<-0;
      }
    }else{
      spatial_grid[[key]]$ca_virus2<-round(site$ca_virus2*exp(-parameters$cavirus_decay/parameters$max_rate));
    }
  }
}

#EJECTION of cell-associated virus to free virus
ca_viral_ejection<-function(i,j,parameters,virus_type=1){
  key<-paste(i,j,sep=",");
  site<-spatial_grid[[key]];
  
  if(virus_type==1 && site$ca_virus1>0){
    num_to_eject<-round(site$ca_virus1*parameters$ca_to_free_proportion);
    if(num_to_eject>0){
      spatial_grid[[key]]$ca_virus1<-site$ca_virus1-num_to_eject;
      spatial_grid[[key]]$free_virus1<-site$free_virus1+num_to_eject;
    }
  }else if(virus_type==2 && site$ca_virus2>0){
    num_to_eject<-round(site$ca_virus2*parameters$ca_to_free_proportion);
    if(num_to_eject>0){
      spatial_grid[[key]]$ca_virus2<-site$ca_virus2-num_to_eject;
      spatial_grid[[key]]$free_virus2<-site$free_virus2+num_to_eject;
    }
  }
}

# OUTPUT 
output_results<-function(settings,results,header){
  output_file<-paste(settings$dirname,'/',settings$out_fname,sep='');
  
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
compute_totals<-function(settings,parameters){
  dynamic$time<<-dynamic$num_updates/parameters$max_rate;
  
  #Count by iterating through spatial_grid
  dynamic$virions1<<-0;
  dynamic$ca_virions1<<-0;
  dynamic$virions2<<-0;
  dynamic$ca_virions2<<-0;
  dynamic$susceptible_cells<<-0;
  dynamic$infected_nonprodcells<<-0;
  dynamic$infected_prodcells<<-0;
  dynamic$dead_cells<<-0;
  dynamic$viral_cells1<<-0;
  dynamic$ca_viral_cells1<<-0;
  dynamic$viral_cells2<<-0;
  dynamic$ca_viral_cells2<<-0;
  dynamic$coinfected_cells<<-0;
  
  for(i in 1:settings$L){
    for(j in 1:settings$L){
      key<-paste(i,j,sep=",");
      site<-spatial_grid[[key]];
      
      dynamic$virions1<<-dynamic$virions1+site$free_virus1;
      dynamic$ca_virions1<<-dynamic$ca_virions1+site$ca_virus1;
      dynamic$virions2<<-dynamic$virions2+site$free_virus2;
      dynamic$ca_virions2<<-dynamic$ca_virions2+site$ca_virus2;
      
      if(site$free_virus1>0) dynamic$viral_cells1<<-dynamic$viral_cells1+1;
      if(site$ca_virus1>0) dynamic$ca_viral_cells1<<-dynamic$ca_viral_cells1+1;
      if(site$free_virus2>0) dynamic$viral_cells2<<-dynamic$viral_cells2+1;
      if(site$ca_virus2>0) dynamic$ca_viral_cells2<<-dynamic$ca_viral_cells2+1;
      
      #Track coinfection: sites with both virus types present (free or CA)
      has_v1<-(site$free_virus1>0 || site$ca_virus1>0);
      has_v2<-(site$free_virus2>0 || site$ca_virus2>0);
      if(has_v1 && has_v2) dynamic$coinfected_cells<<-dynamic$coinfected_cells+1;
      
      if(site$cell_state==1) dynamic$susceptible_cells<<-dynamic$susceptible_cells+1;
      if(site$cell_state==2) dynamic$infected_nonprodcells<<-dynamic$infected_nonprodcells+1;
      if(site$cell_state==3) dynamic$infected_prodcells<<-dynamic$infected_prodcells+1;
      if(site$cell_state==4) dynamic$dead_cells<<-dynamic$dead_cells+1;
    }
  }
  
  total_virions<-dynamic$virions1+dynamic$virions2;
  if (total_virions > dynamic$max_vl) {
    dynamic$max_vl <<-  total_virions;
    dynamic$max_vl_time <<-  dynamic$time;
  }
}

#END of simulation - save grids
end_simulation<-function(settings,parameters){
  if(settings$save_grids){
    #Extract arrays from spatial_grid structure for saving
    cell_state_array<-matrix(0,nrow=settings$L,ncol=settings$L);
    free_virus1_array<-matrix(0,nrow=settings$L,ncol=settings$L);
    ca_virus1_array<-matrix(0,nrow=settings$L,ncol=settings$L);
    free_virus2_array<-matrix(0,nrow=settings$L,ncol=settings$L);
    ca_virus2_array<-matrix(0,nrow=settings$L,ncol=settings$L);
    
    for(i in 1:settings$L){
      for(j in 1:settings$L){
        key<-paste(i,j,sep=",");
        cell_state_array[i,j]<-spatial_grid[[key]]$cell_state;
        free_virus1_array[i,j]<-spatial_grid[[key]]$free_virus1;
        ca_virus1_array[i,j]<-spatial_grid[[key]]$ca_virus1;
        free_virus2_array[i,j]<-spatial_grid[[key]]$free_virus2;
        ca_virus2_array[i,j]<-spatial_grid[[key]]$ca_virus2;
      }
    }
    save(cell_state_array,file=paste(settings$dirname,'/cell_state',sep=''));
    save(free_virus1_array,file=paste(settings$dirname,'/free_virus1',sep=''));
    save(ca_virus1_array,file=paste(settings$dirname,'/ca_virus1',sep=''));
    save(free_virus2_array,file=paste(settings$dirname,'/free_virus2',sep=''));
    save(ca_virus2_array,file=paste(settings$dirname,'/ca_virus2',sep=''));
  }
  if(settings$visualize=='pdf'){
    dev.off();
  }
  
  #Plot dynamics
  plot_dynamics(settings);
  
  print(paste("Highest log VL",ifelse(dynamic$max_vl>0,log10(dynamic$max_vl),0),"at t=",dynamic$max_vl_time));
  print(paste("End of simulation at t=",dynamic$time));
}