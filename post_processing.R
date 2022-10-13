## USDOS Results Processing ##
# Generates plots, maps, summaries, and data files from USDOS summary and detail output files

# Code by Deedra Murrieta & Katie Owers, edits by Lauren Smith & Sophie McKee

# Summary and Detail files must be saved in the directory "Files_To_Process" in the Post_Processing directory.


#####################################################################################################

#' Process USDOS output files
#'
#' Generate summaries of USDOS runs as figures, maps, files, and tables based on the model's summary and detail file outputs. Set the working directory as a folder containing the output files to be processed. 
#'
#' @param  export.datafiles Export the files used to calculate USDOS metrics? 0 = no (default), 1 = export wide files that contain each run's values in a column, followed by a column with the run type, 2 = export long files hat contain all run values in a single column wit another column of tun types, 3 = export both wide and long files
#' @param  results.report Should an html results report explaining results and plots be generated? (TRUE/FALSE, default = TRUE)
#' @param  summaryTable Generate a summary table for Duration, animals infected, premises infected, and epidemic extent? (TRUE/FALSE, default = TRUE)
#' @param  duration Generate figures for the duration of infection metric (as measured in model timesteps--days for FMD). (TRUE/FALSE, default = TRUE)
#' @param  Dur_min A mimumum duration of interest (default = 13)
#' @param  Dur_cutoff A duration above which an outreak would be considered "large" (default = 100)
#' @param  premInf Generate figures for the infected premises metric? (TRUE/FALSE, default = TRUE)
#' @param  PremInf_min A minimum number of premises infected of interest (default = 10)
#' @param  PremInf_cutoff A number of infected premises above which an outbreak would be considered 'large' (default = 1000)
#' @param  premReport = Generate figures for the reported premises metric? 
#' @param  ReportedPrems_min A minimum number of premises reported that is of interest (default = 5)
#' @param  ReportedPrems_cutoff A number of infected premises above which an outbreak would be considered 'large' (default = 100)
#' @param  epidemicExtent Generate figures for the epidemic extent metric? (TRUE/FALSE, default = TRUE)
#' @param  EpidExt_min A minimum number of counties infected that is of interest (default = 1)
#' @param  EpidExt_cutoff A number of infected counties above which an outbreak would be considered 'large'
#' @param  movementBan  Generate figures for the number of geographies (states or counties as specified in the configuration file)? (TRUE/FALSE, default = TRUE)
#' @param  MB_min A minimum number of areas affected by a movement ban that would be of interest (default= 0)
#' @param  MB_cutoff A number of areas affected by a movement ban above which an outbreak would be considered 'large'
#' @param  premisesCulled Generate figures for the premises culled metric? Only used for runs with a control type that includes culling. (TRUE/FALSE, default = TRUE)
#' @param  PremCull_min A minimumm number of premises culled that is of interest (default = 0)
#' @param  PremCull_cutoff A number of culled premises above which an outbreak would be considered 'large' (default = 10)
#' @param  premisesVax Generate figures for the premises vaccinated metric? Only used for runs with a control type that includes vaccination. (TRUE/FALSE, default = TRUE)
#' @param  PremVax_min A minimum number of premises vaccinated that is of interest (default = 0)
#' @param  PremVax_cutoff A number of premises vaccinated above which an outbreak would be consdiered 'large'
#' @param  animalsInfected Genreate figures for the number of animals infected metric? (TRUE/FALSE, default = TRUE)
#' @param  Anim_min A minimum number of animals infected that is of interest (default = 10)
#' @param  Anim_cutoff A number of animals infected above which an outbreak would be considered 'large' (default = 1000)
#' @param  countyRisk Generate figures for the county risk metric? (TRUE/FALSE, default = TRUE)
#' @param  CountyRisk_min A minimum county risk that is of interest (default = 0)
#' @param  CountyRisk_cutoff A county risk that is considered large (default = 0.001, crresponding to a county becoming infected when outbreaks are seeded in at least four other counties. (TRUE/FALSE, default = TRUE)
#' @param  localSpread Generate figures for the proportion of spread that is due to local spread (versus shipment). (TRUE/FALSE, default = TRUE)
#' @param  color_palette The color palette for maps other than local spread (see ls_match option below to change the local spread map). By default this is "color_red" (the RColorBrewer 'OrRd' palette). Other options are "color_blue" (RColorBrewer 'Blues'), "color_orange" ('Oranges'), "color_pink" ('RdPu'), "color_bluepurple" ('BuPu'), and "color_yellow" (custom)
#' @param  ls_match Should the local spread map's color palette match the other maps (TRUE, the default value) or use its default of color_bluepurple (FALSE)

#' @details 
#' [metric]_min Allows users to plot results excluding simulations below a certain minimum value. For example, since the duration of outbreaks that do not spread is 13 days, setting a Dur_minimal of 13 would generate plots which exclude outbreaks that do not persist beyond the index infection.
#' [metric]_cutoff Allows users to set a cutoff over which an outbreak would be considered "large". Outbreaks that go above this threshold will be plotted separately fromthose that stay below it


#' @author Send bug reports, suggestions, corrections, or comments to webblaboratory(at)gmail.com
#' 
#####################################################################################################
#####################################################################################################

#####     Function       #######

#####################################################################################################
#####################################################################################################

processUSDOS = function(export.datafiles= 0,
                        results.report = FALSE,
                        summaryTable = TRUE,
                        duration = TRUE,
                        Dur_min = 13,
                        Dur_cutoff = 100,
                        premInf = TRUE,
                        PremInf_min = 10, 
                        PremInf_cutoff = 1000,
                        premReport = TRUE,
                        ReportedPrems_min = 5,
                        ReportedPrems_cutoff = 100,
                        epidemicExtent = TRUE,
                        EpidExt_min= 1,
                        EpidExt_cutoff=50,
                        movementBan = TRUE,
                        MB_min = 0,
                        MB_cutoff = 1,
                        premisesCulled = TRUE,
                        PremCull_min = 0,
                        PremCull_cutoff = 10,
                        premisesVax = TRUE,
                        PremVax_min = 0,
                        PremVax_cutoff = 2,
                        animalsInfected = TRUE,
                        Anim_min = 10,
                        Anim_cutoff = 1000,
                        countyRisk = TRUE,
                        CountyRisk_min = 0,
                        CountyRisk_cutoff = 0.001,
                        localSpread = TRUE,
                        color_palette = "color_red",
                        ls_match = TRUE
                        
)
{
  
  
  #####################################################################################################
  #####################################################################################################
  
  #####     Setup      #######
  
  #####################################################################################################
  #####################################################################################################
  
  ## Load Packages ##
  library(stats); library(maps);  library(mapdata); library(tidyr); library(fields); library(foreach);
  library(dplyr); library(RColorBrewer); library(rgdal); library(reshape2); library(reshape); library(data.table); library(knitr);
  #library(ggpubr); 
  library(ggplot2); library(ggmap); library(kableExtra); library(magrittr)
  library(tidyverse)
  
  # This function finds the directory in which the the current file is located.
  getCurrentFileLocation <-  function()
  {
    this_file <- commandArgs() %>% 
      tibble::enframe(name = NULL) %>%
      tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
      dplyr::filter(key == "--file") %>%
      dplyr::pull(value)
    if (length(this_file)==0)
    {
      this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
  }
  
  path0 <- getCurrentFileLocation()
  setwd(path0)
  
  ## Source map_by_fips 
  source("map_by_fips_standalone.R")
  
  ## Generating color palettes and generic scales ##
  cbPalette <- c("#D55E00", "#CC79A7", "#56B4E9", "#009E73", "#0072B2", "#000000", "#F0E442", "#E69F00") # Violin plots
  
  if (color_palette == "color_red") {palette = brewer.pal(8, "OrRd")
  } else if(color_palette == "color_bluepurple"){palette = brewer.pal(8, "BuPu")
  } else if(color_palette == "color_blue"){palette = brewer.pal(8, "Blues")
  } else if(color_palette == "color_yellow"){palette = c("ivory3","lightgoldenrodyellow","lightgoldenrod", "gold2", "goldenrod", "darkgoldenrod") 
  } else if(color_palette == "color_pink"){palette = brewer.pal(8, "RdPu")
  } else if(color_palette == "color_orange"){palette = brewer.pal(8, "Oranges")
  } else if(color_palette == "color_YlOrRd"){palette = rev(c("#000000", "#700000", "#8F0000", "#B10000", "#CE2900", "#E25700", "#F37B00", "#FF9C08", "#FFB954", "#FFD47F", "#FFEBA3", "#FFFBBF"))}
  
  
  ## Set up geography files ##
  county.summary <- county.fips
  county.summary[which(county.summary$fips==46113),1]<- 46102 # South Dakota county that changed name/number 
  county.summary[nrow(county.summary)+1,] <- c(51550, "virginia,chesapeake") # Adding in a county missing from county.fips
  county.summary$fips<-as.numeric(county.summary$fips)
  county.summary <- separate(county.summary, polyname, into = c("county", "extra"), sep = ":", fill = "right")
  county.summary$county[county.summary$county == 'montana,park'] <- "montana,yellowstone national"
  county.summary <- county.summary %>%
    select(fips, county) %>%
    dplyr::rename('polyname' = county) %>%
    distinct()
  
  county.summary$polyname=as.factor(county.summary$polyname)
  
  
  #### Read and format Summary Files ####
  
  ## Identify all the summary files ##
  pathfiles <- file.path(path0, "Files_To_Process/")
  
  summary.files <- list.files(path = pathfiles, recursive = TRUE, pattern = "_summary.txt", full.names = FALSE)
  
  ## Read in all the summary files ##
  ### Storing as a single file with many columns    
  
  # Initiate A vector to fill with the run types identified from file names
  run.types=NULL 
  
  # Summary file reading
  ### This piece takes several minutes
  
  
  for(i in 1:length(summary.files)){
    summary.file <- summary.files[i]
    summary.file.path <- paste(pathfiles, summary.file, sep = "")
    summary.res <- fread(summary.file.path, header = TRUE)
    summary.res <- as.data.frame(summary.res)
    summary.res2 <- summary.res[,-grep( "RunTime" , colnames( summary.res ))] 
    summary.res2$Type <- unlist(strsplit(summary.file, "_FLAPS"))[1]
    run.types=c(run.types,summary.res2$Type[1])
    names(summary.res2)[names(summary.res2) == 'Seed_FIPS'] <- "fips"
    county.summary <- merge(county.summary, summary.res2, by = "fips", all = TRUE)
  }
  
  
  rm(summary.res,summary.res2)
  
  # Identify run types included
  run.types=unique(run.types);run.types
  
  ## To check for <100 runs per type and account for possibility that users run >100 iterations per county, calculate how many they ran. 
  runs_per_ctrl_type= length(summary.files)/length(run.types)
  
  # Error if <100 summary files per run type
  if( runs_per_ctrl_type<100 ) c(warning('At least 100 iterations are required to capture uncertainty in model predictions. Do not analyze fewer than 100 runs.'), rm(county.summary))
  
  
  # Error if summary files per run type is not a whole number
  if( runs_per_ctrl_type!=round(runs_per_ctrl_type) ) c(warning('Number of summary files per run type is not a whole number. Check data files for missing or extra data files.'), 
                                                        rm(county.summary))
  
  
 
  
  #####################################################################################################
  #####################################################################################################
  
  ##### Duration  (Summary Files)      #######
  
  #####################################################################################################
  #####################################################################################################
  
  path_output <- file.path(path0, "Output_Files/")
  dir.create(path_output, showWarnings = FALSE)
  
  
  if (duration == TRUE){
    setwd(path_output)
    Dur=county.summary[,grepl( "fips|polyname|Duration|Type" , names( county.summary ) )]
    names(Dur)[seq(4,ncol(Dur),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(Dur))/2-1))))
    names(Dur)[seq(3,ncol(Dur),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(Dur))/2-1)))) 
    
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(Dur, "Duration.csv",row.names=F)} 

    
    
    # Convert to long shape for plotting
    Dur.long <- reshape(Dur, idvar = c("fips", "polyname"), direction = "long",
                        v.names = c("type", "Value"),
                        varying = list(c(grep("type_", colnames(Dur))),
                                       c(grep("run_", colnames(Dur)))))
    
    rownames(Dur.long) <- seq(1:nrow(Dur.long))
    Dur.long <- Dur.long[,-3]
    Dur.long <- Dur.long[!is.na(Dur.long$Value),]
    
    if (export.datafiles == 2 | export.datafiles == 3) {write.csv(Dur.long, "Duration_long.csv",row.names=F)}
    
    ## Divide the data into subsets for plotting/mapping ##
    
    Dur.long.overMin <- Dur.long[Dur.long$Value > Dur_min,]
    Dur.long.low <- Dur.long[Dur.long$Value <= Dur_cutoff,]
    Dur.long.high <- Dur.long[Dur.long$Value > Dur_cutoff,]
    
    
    ## Histograms ##
    jpeg("Duration_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(Dur.long.high$Value, breaks = 25, xlab = "Duration", 
         main = paste0("Duration >", as.character(Dur_cutoff)," timesteps"))
    dev.off()
    
    jpeg("Duration_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(Dur.long.low$Value, breaks = 25, xlab = "Duration", 
         main = paste0("Duration \u2264", as.character(Dur_cutoff)," timesteps"))
    dev.off()
    
    jpeg("Duration_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(Dur.long.overMin$Value, breaks = 25, xlab = "Duration", 
         main = paste0("Duration >", as.character(Dur_min)," timesteps"))
    dev.off()
    
    ## violin plots ##
    jpeg("Duration_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = Dur.long.overMin, aes(x = Dur.long.overMin$type, y = Dur.long.overMin$Value, color = Dur.long.overMin$type)) + 
        scale_x_discrete(limits = levels(Dur.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "Duration of outbreak", 
             title = eval(substitute(paste0("Duration >",v," timesteps"), list(v=Dur_min)))) + 
        scale_color_manual(values = cbPalette) 
    })
    dev.off()
    
    jpeg("Duration_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = Dur.long.low, aes(x = Dur.long.low$type, y = Dur.long.low$Value, color = Dur.long.low$type)) + 
        scale_x_discrete(limits = levels(Dur.long.low$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        scale_color_manual(values = cbPalette) + 
        labs(x = "Run Type", y = "Duration of outbreak",
             title=eval(substitute(paste0("Duration \u2264",v," timesteps"), list(v=Dur_cutoff)))) 
    })
    dev.off()
    
    jpeg("Duration_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = Dur.long.high, aes(x = Dur.long.high$type, y = Dur.long.high$Value, color = Dur.long.high$type)) + 
        scale_x_discrete(limits = levels(Dur.long.high$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +  
        scale_color_manual(values = cbPalette) + 
        labs(x = "Run Type", y = "Duration of outbreak",
             title=eval(substitute(paste0("Duration >",v, " timesteps"), list(v=Dur_cutoff))))
    })
    dev.off()
    
    
    
    Dur_Median_RepType=tapply(Dur.long$Value,Dur.long$type,FUN=median)
    Dur_Upper_RepType<-tapply(Dur.long$Value,Dur.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
    
    ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
    # two-column dataframes with fips for mapping
    
    Dur_medians=NULL
    Dur_uppers=NULL
    
    for (i in 1:length(run.types)){
      df <- Dur[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
      df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
      newname <- df1[1,4]
      df.new <- df[,-grep("type_", colnames(df))]
      
      # Calculate median and upper 2.5%
      df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
      df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
      
      # Create the vectors for map scales
      Dur_medians=c(Dur_medians,df.new$median)
      Dur_uppers=c(Dur_uppers,df.new$upper)
      
      # save datasets for mapping, one for each level with medians and uppers with name of control type
      assign(paste0("Dur_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
      assign(paste0("Dur_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
      
      rm(df, df1,df.new,newname)
    }
    
    ## Median maps
    
    # Create a vector for scale for median maps
    Dur_median_values <- round(c(0, Dur_min, mean(Dur_medians, na.rm = TRUE), 
                                 mean(Dur_medians, na.rm = TRUE) + sd(Dur_medians, na.rm = TRUE), 
                                 mean(Dur_medians, na.rm = TRUE) + 3*sd(Dur_medians, na.rm = TRUE), 
                                 mean(Dur_medians, na.rm = TRUE) + 6*sd(Dur_medians, na.rm = TRUE), 
                                 max(Dur_medians, na.rm = TRUE)), 0)
    Dur_median_values<-unique(sort(Dur_median_values));Dur_median_values
    if (length(Dur_median_values) ==2){Dur_median_values <- c(0, max(Dur_median_values)-1, max(Dur_median_values) + 1)}
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("Dur_med_",run.types[i])) 
      name_df <- na.omit(name_df)
      jpeg(paste0(path_output, paste0("Duration_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = Dur_median_values, color.sequence = palette, 
                  legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
      dev.off()
    }
    
    ## Upper 2.5% maps
    
    # Create a vector for scale for upper maps
    Dur_upper_values <- round(c(0, Dur_min, mean(Dur_uppers, na.rm = TRUE), 
                                mean(Dur_uppers, na.rm = TRUE) + sd(Dur_uppers, na.rm = TRUE), 
                                mean(Dur_uppers, na.rm = TRUE) + 3*sd(Dur_uppers, na.rm = TRUE), 
                                mean(Dur_uppers, na.rm = TRUE) + 6*sd(Dur_uppers, na.rm = TRUE), 
                                max(Dur_uppers, na.rm = TRUE)), 0)
    if(6*sd(Dur_uppers, na.rm = TRUE)> max(Dur_uppers, na.rm = TRUE)) Dur_upper_values<-Dur_upper_values[c(1:5,7)]
    Dur_upper_values<-unique(sort(Dur_upper_values));Dur_upper_values
    
    
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("Dur_upper_",run.types[i])) 
      jpeg(paste0(path_output, paste0("Duration_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = Dur_upper_values, color.sequence = palette, 
                  legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
      dev.off()
    }
  }
  ##################################################################
  ##################################################################
  
  #### Number of Premises Infected (Summary Files) ####
  
  ##################################################################
  ##################################################################
  
  if (premInf == TRUE){
    PremInf=county.summary[,grepl( "fips|polyname|Num_Inf|Type" , names( county.summary ) )]
    names(PremInf)[seq(4,ncol(PremInf),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(PremInf))/2-1))))
    names(PremInf)[seq(3,ncol(PremInf),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(PremInf))/2-1)))) 
    
    setwd(path_output)
    
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(PremInf, "PremisesInfected.csv",row.names=F)}
    
    # Convert to long format for plotting
    PremInf.long <- reshape(PremInf, idvar = c("fips", "polyname"), direction = "long", 
                            v.names = c("type", "Value"), 
                            varying = list(c(grep("type_", colnames(PremInf))), 
                                           c(grep("run_", colnames(PremInf)))))
    
    rownames(PremInf.long) <- seq(1:nrow(PremInf.long))
    PremInf.long <- PremInf.long[,-3]
    PremInf.long <- PremInf.long[!is.na(PremInf.long$Value),]
    
    if (export.datafiles == 2 | export.datafiles == 3) {write.csv(PremInf.long,"PremisesInfected_long.csv",row.names = F)}

    
    sum(PremInf.long$Value>1)/length(PremInf.long$Value)
    
    ## Divide the data into subsets for plotting ## 
    
    PremInf.long.overMin <- PremInf.long[PremInf.long$Value > PremInf_min,]
    PremInf.long.low <- PremInf.long[PremInf.long$Value <= PremInf_cutoff,]
    PremInf.long.high <- PremInf.long[PremInf.long$Value > PremInf_cutoff,]
    
    ## Histograms ##
    jpeg("PremInf_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(PremInf.long.high$Value, breaks = 100, xlab = "# Infected Premises",
         main = paste0("# Infected Premises >", as.character(PremInf_cutoff)))
    dev.off()
    
    jpeg("PremInf_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(PremInf.long.low$Value, breaks = 100,xlab = "# Infected Premises",
         main = paste0("# Infected Premises \u2264", as.character(PremInf_cutoff)))
    dev.off()
 
    jpeg("PremInf_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(PremInf.long.overMin$Value, breaks = 100, 
         xlab = "# Infected Premises", main = paste0("# Infected Premises >", as.character(PremInf_min)))
    dev.off()
    
    ## Violin Plots ##
    jpeg("PremInf_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = PremInf.long.overMin, aes(x = PremInf.long.overMin$type, y = PremInf.long.overMin$Value, color = PremInf.long.overMin$type)) + 
        scale_x_discrete(limits = levels(PremInf.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "# Infected Premises", 
             title = eval(substitute(paste0("# Infected premises >",v), list(v=PremInf_min)))) +  
        scale_color_manual(values = cbPalette)
    })
    dev.off()
    
    jpeg("PremInf_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = PremInf.long.low, aes(x = PremInf.long.low$type, y = PremInf.long.low$Value, color = PremInf.long.low$type)) + 
        scale_x_discrete(limits = levels(PremInf.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "# Infected Premises", title = paste0("# Infected premises \u2264", as.character(PremInf_cutoff) )) + 
        scale_color_manual(values = cbPalette)
    })
    dev.off()
    
    jpeg("PremInf_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = PremInf.long.high, aes(x = PremInf.long.high$type, y = PremInf.long.high$Value, color = PremInf.long.high$type)) + 
        scale_x_discrete(limits = levels(PremInf.long.overMin$type)) +  theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "# Infected Premises", 
             title = eval(substitute(paste0("# Infected premises >",v), list(v=PremInf_cutoff)))) + 
        scale_color_manual(values = cbPalette)
    })
    dev.off()
    
    
    
    PremInf_Median_RepType=tapply(PremInf.long$Value,PremInf.long$type,FUN=median)
    PremInf_Upper_RepType<-tapply(PremInf.long$Value,PremInf.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
    
    ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
    # two-column dataframes with fips for mapping
    
    PremInf_medians=NULL
    PremInf_uppers=NULL
    gr10IPs=NULL
    
    for (i in 1:length(run.types)){
      
      df <- PremInf[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
      df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
      newname <- df1[1,4]
      df.new <- df[,-grep("type_", colnames(df))]
      
      # Calculate median and upper 2.5%
      df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
      df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
      
      # Calculate proportion of runs with >10 IPs
      df.new$gr10 = apply(df.new[grep("run_", colnames(df.new))], 1, function(x) length(x[x>10]))
      gr10IPs[i]=sum(df.new$gr10, na.rm = TRUE)
      
      # Create the vectors for map scales
      PremInf_medians=c(PremInf_medians,df.new$median)
      PremInf_uppers=c(PremInf_uppers,df.new$upper)
      
      # save datasets for mapping, one for each level with medians and uppers with name of control type
      assign(paste0("PremInf_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
      assign(paste0("PremInf_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
      
      rm(df, df1,df.new,newname)
    }
    
    # For Summary Table - find the total number of times it was >10 across the counties, then divide by total possible (304900) and * 100 for percentage
    gr10IPs <- gr10IPs/3049
    
    ## Median maps
    
    # Create a vector for scale for median maps
    PremInf_median_values <- round(c(0, 1, PremInf_min, mean(PremInf_medians, na.rm = TRUE), 
                                     mean(PremInf_medians, na.rm = TRUE) + sd(PremInf_medians, na.rm = TRUE), 
                                     mean(PremInf_medians, na.rm = TRUE) + 3*sd(PremInf_medians, na.rm = TRUE), 
                                     mean(PremInf_medians, na.rm = TRUE) + 6*sd(PremInf_medians, na.rm = TRUE), 
                                     max(PremInf_medians, na.rm = TRUE)), 0)
    PremInf_median_values<-unique(sort(PremInf_median_values))
    
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("PremInf_med_",run.types[i])) 
      if (var(name_df$median, na.rm=TRUE)==0) {name_df$median[1]<-(name_df$median[1]+0.001)} # map_by_fips doesn't work with all zeros, so add a tiny bit.
      jpeg(paste0(path_output, paste0("PremInf_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = PremInf_median_values, color.sequence = palette, 
                  legend.spacing = 4, legend.shrink = 0.5, legend.width = 1)
      dev.off()
    }
    
    ## Upper 2.5% maps
    
    # Create a vector for scale for upper maps
    PremInf_upper_values <- round(c(0, 1, PremInf_min, mean(PremInf_uppers, na.rm = TRUE), 
                                    mean(PremInf_uppers, na.rm = TRUE) + sd(PremInf_uppers, na.rm = TRUE), 
                                    mean(PremInf_uppers, na.rm = TRUE) + 3*sd(PremInf_uppers, na.rm = TRUE), 
                                    mean(PremInf_uppers, na.rm = TRUE) + 6*sd(PremInf_uppers, na.rm = TRUE), 
                                    max(PremInf_uppers, na.rm = TRUE)), 0)
    if(6*sd(PremInf_uppers, na.rm = TRUE)> max(PremInf_uppers, na.rm = TRUE)) PremInf_upper_values<-PremInf_upper_values[c(1:5,7)]
    PremInf_upper_values<-unique(sort(PremInf_upper_values))
    
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("PremInf_upper_",run.types[i])) 
      jpeg(paste0(path_output, paste0("PremInf_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = PremInf_upper_values, color.sequence = palette, 
                  legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
      dev.off()
    }
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #### Number of Reported Premises (Summary Files) ####
  
  #####################################################################################################
  #####################################################################################################
  
  if (premReport == TRUE){
    
    ReportedPrems <- county.summary[,grepl( "fips|polyname|Num_Reports|Type" , names( county.summary ) )]
    
    deleteme=c(0,rep(NA,ncol(ReportedPrems)-1))
    
    # this gets rid of adjacent "type" columns, which indicate that a movement ban wasn't implemented
    for(i in 2:ncol(ReportedPrems)){
      deleteme[i]=ifelse(grepl("Type", colnames(ReportedPrems)[i-1]) & grepl("Type", colnames(ReportedPrems)[i])|
                           grepl("Type", colnames(ReportedPrems)[i]) & grepl("Type", colnames(ReportedPrems)[i+1]),1,0)
    }
    
    ReportedPrems=ReportedPrems[,deleteme==0]
    
    if (ncol(ReportedPrems)>2){
      setwd(path_output)
      
      names(ReportedPrems)[seq(4,ncol(ReportedPrems),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(ReportedPrems))/2-1))))
      names(ReportedPrems)[seq(3,ncol(ReportedPrems),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(ReportedPrems))/2-1))))
      
      if (export.datafiles == 1 | export.datafiles == 3) {write.csv(ReportedPrems, "ReportedPremises.csv",row.names=F)}

      # Convert to long shape for plotting  
      ReportedPrems.long <- reshape(ReportedPrems, idvar = c("fips", "polyname"), direction = "long", 
                                    v.names = c("type", "Value"), 
                                    varying = list(c(grep("type_", colnames(ReportedPrems))), 
                                                   c(grep("run_", colnames(ReportedPrems)))))
      
      rownames(ReportedPrems.long) <- seq(1:nrow(ReportedPrems.long))
      ReportedPrems.long <- ReportedPrems.long[,-3]
      ReportedPrems.long <- ReportedPrems.long[!is.na(ReportedPrems.long$Value),] 
      if (export.datafiles == 2 | export.datafiles == 3) {write.csv(ReportedPrems.long,"ReportedPremises_long.csv",row.names = F)}
      
      ReportedPrems.long.overMin <- ReportedPrems.long[ReportedPrems.long$Value > ReportedPrems_min,]
      ReportedPrems.long.low <- ReportedPrems.long[ReportedPrems.long$Value <= ReportedPrems_cutoff,]
      ReportedPrems.long.high <- ReportedPrems.long[ReportedPrems.long$Value > ReportedPrems_cutoff,]
      
      ## Histograms ##
      jpeg("ReportedPrems_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(ReportedPrems.long.high$Value, breaks = 100, xlab = "# Reported Premises",
           main = paste0("# Reported Premises >", as.character(ReportedPrems_cutoff)))
      dev.off()
      
      jpeg("ReportedPrems_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(ReportedPrems.long.low$Value, breaks = 100,xlab = "# Reported Premises",
           main = paste0("# Reported Premises \u2264", as.character(ReportedPrems_cutoff)))
      dev.off()
      
      jpeg("ReportedPrems_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(ReportedPrems.long.overMin$Value, breaks = 100, 
           xlab = "# Reported Premises", main = paste0("# Reported Premises >", as.character(ReportedPrems_min)))
      dev.off()
      
      ## Violin Plots ##
      jpeg("ReportedPrems_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = ReportedPrems.long.overMin, aes(x = ReportedPrems.long.overMin$type, y = ReportedPrems.long.overMin$Value, color = ReportedPrems.long.overMin$type)) + 
          scale_x_discrete(limits = levels(ReportedPrems.long.overMin$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Reported Premises", 
               title = eval(substitute(paste0("# Reported premises >",v), list(v=ReportedPrems_min)))) +  
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      jpeg("ReportedPrems_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = ReportedPrems.long.low, aes(x = ReportedPrems.long.low$type, y = ReportedPrems.long.low$Value, color = ReportedPrems.long.low$type)) + 
          scale_x_discrete(limits = levels(ReportedPrems.long.overMin$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Reported Premises", title = paste0("# Reported premises \u2264", as.character(ReportedPrems_cutoff) )) + 
          scale_color_manual(values = cbPalette)
        dev.off()
        
        jpeg("ReportedPrems_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
        ggplot() + geom_violin(data = ReportedPrems.long.high, aes(x = ReportedPrems.long.high$type, y = ReportedPrems.long.high$Value, color = ReportedPrems.long.high$type)) + 
          scale_x_discrete(limits = levels(ReportedPrems.long.overMin$type)) +  theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Reported Premises", 
               title = eval(substitute(paste0("# Reported premises >",v), list(v=ReportedPrems_cutoff)))) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      
      ReportedPrems_Median_RepType=tapply(ReportedPrems.long$Value,ReportedPrems.long$type,FUN=median)
      ReportedPrems_Upper_RepType<-tapply(ReportedPrems.long$Value,ReportedPrems.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
      
      
      ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
      # two-column dataframes with fips for mapping
      
      ReportedPrems_medians=NULL
      ReportedPrems_uppers=NULL
      
      for (i in 1:length(run.types)){
        
        df <- ReportedPrems[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
        df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
        newname <- df1[1,4]
        df.new <- df[,-grep("type_", colnames(df))]
        
        # Calculate median and upper 2.5%
        df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
        df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
        
        # Create the vectors for map scales
        ReportedPrems_medians=c(ReportedPrems_medians,df.new$median)
        ReportedPrems_uppers=c(ReportedPrems_uppers,df.new$upper)
        
        # save datasets for mapping, one for each level with medians and uppers with name of control type
        assign(paste0("ReportedPrems_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
        assign(paste0("ReportedPrems_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
        
        rm(df, df1,df.new,newname)
      }
      
      ## Median maps
      
      # Create a vector for scale for median maps
      ReportedPrems_median_values <- round(c(0, 1, ReportedPrems_min, mean(ReportedPrems_medians, na.rm = TRUE), 
                                             mean(ReportedPrems_medians, na.rm = TRUE) + sd(ReportedPrems_medians, na.rm = TRUE), 
                                             mean(ReportedPrems_medians, na.rm = TRUE) + 3*sd(ReportedPrems_medians, na.rm = TRUE), 
                                             mean(ReportedPrems_medians, na.rm = TRUE) + 6*sd(ReportedPrems_medians, na.rm = TRUE), 
                                             max(ReportedPrems_medians, na.rm = TRUE)), 0)
      ReportedPrems_median_values<-unique(sort(ReportedPrems_median_values))
      
      
      for (i in 1:length(run.types)){
        name_df=get(paste0("ReportedPrems_med_",run.types[i])) 
        jpeg(paste0(path_output, paste0("ReportedPrems_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = ReportedPrems_median_values, color.sequence = palette, 
                    legend.spacing = 4, legend.shrink = 0.5, legend.width = 1)
        dev.off()
      }
      
      ## Upper 2.5% maps
      
      # Create a vector for scale for upper maps
      ReportedPrems_upper_values <- round(c(0, 1, ReportedPrems_min, mean(ReportedPrems_uppers, na.rm = TRUE), 
                                            mean(ReportedPrems_uppers, na.rm = TRUE) + sd(ReportedPrems_uppers, na.rm = TRUE), 
                                            mean(ReportedPrems_uppers, na.rm = TRUE) + 3*sd(ReportedPrems_uppers, na.rm = TRUE), 
                                            mean(ReportedPrems_uppers, na.rm = TRUE) + 6*sd(ReportedPrems_uppers, na.rm = TRUE), 
                                            max(ReportedPrems_uppers, na.rm = TRUE)), 0)
      if(6*sd(ReportedPrems_uppers, na.rm = TRUE)> max(ReportedPrems_uppers, na.rm = TRUE)) ReportedPrems_upper_values<-ReportedPrems_upper_values[c(1:5,7)]
      ReportedPrems_upper_values<-unique(sort(ReportedPrems_upper_values))
      
      
      for (i in 1:length(run.types)){
        name_df=get(paste0("ReportedPrems_upper_",run.types[i])) 
        jpeg(paste0(path_output, paste0("ReportedPrems_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = ReportedPrems_upper_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
        dev.off()
      }
      
    } 
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #### Number of affected counties aka Epidemic Extent (nAffCounties in Summary Files) ####
  
  #####################################################################################################
  #####################################################################################################
  
  if (epidemicExtent == TRUE){
    EpidExt=county.summary[,grepl( "fips|polyname|nAffCounties|Type" , names( county.summary ) )]
    names(EpidExt)[seq(4,ncol(EpidExt),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(EpidExt))/2-1))))
    names(EpidExt)[seq(3,ncol(EpidExt),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(EpidExt))/2-1)))) 
    
    setwd(path_output)
    
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(EpidExt, "EpidemicExtent.csv",row.names=F)}
  
    
    # Convert to long shape for plotting  
    EpidExt.long <- reshape(EpidExt, idvar = c("fips", "polyname"), direction = "long", 
                            v.names = c("type", "Value"), 
                            varying = list(c(grep("type_", colnames(EpidExt))), 
                                           c(grep("run_", colnames(EpidExt)))))
    
    rownames(EpidExt.long) <- seq(1:nrow(EpidExt.long))
    EpidExt.long <- EpidExt.long[,-3]
    EpidExt.long <- EpidExt.long[!is.na(EpidExt.long$Value),] 
    if (export.datafiles == 2 | export.datafiles == 3) {write.csv(EpidExt.long,"EpidemicExtent_long.csv",row.names = F)}

    
    summary(EpidExt.long$Value)
    
    ## Divide the data into subsets for plotting ## 
    
    EpidExt.long.overMin <- EpidExt.long[EpidExt.long$Value> EpidExt_min,]
    EpidExt.long.low <- EpidExt.long[EpidExt.long$Value <= EpidExt_cutoff,]
    EpidExt.long.high <- EpidExt.long[EpidExt.long$Value > EpidExt_cutoff,]
    
    ## Histograms ##
    jpeg("EpidExt_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(EpidExt.long.high$Value, breaks = 25, xlab = "# Infected Counties", 
         main = paste0("# Infected counties > ",as.character(EpidExt_cutoff)))
    dev.off()
    
    jpeg("EpidExt_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(EpidExt.long$Value[ EpidExt.long$Value <EpidExt_cutoff], breaks = 25, 
         xlab = "# Infected Counties", 
         main = paste0("# Infected counties \u2264",as.character(EpidExt_cutoff)))
    dev.off()
    

    jpeg("EpidExt_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(EpidExt.long$Value[EpidExt.long$Value > 1], breaks = 25, xlab = "# Infected Counties", 
         main = paste0("# Infected counties > ",as.character(EpidExt_min)))
    dev.off()
    
    ## Violin plots ##
    jpeg("EpidExt_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = EpidExt.long.overMin, aes(x = EpidExt.long.overMin$type, y = EpidExt.long.overMin$Value, color = EpidExt.long.overMin$type)) + 
        scale_x_discrete(limits = levels(EpidExt.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "# Infected Premises", 
             title = eval(substitute(paste0("# Infected counties >",v), list(v=EpidExt_min)))) +  
        scale_color_manual(values = cbPalette) 
    })
    dev.off()
    
    jpeg("EpidExt_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = EpidExt.long.low, aes(x = EpidExt.long.low$type, y = EpidExt.long.low$Value, color = EpidExt.long.low$type)) + 
        scale_x_discrete(limits = levels(EpidExt.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "# Infected Counties", title = paste0("# Infected counties \u2264", as.character(EpidExt_cutoff) )) + 
        scale_color_manual(values = cbPalette)
    })
    dev.off()
    
    jpeg("EpidExt_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = EpidExt.long.high, aes(x = EpidExt.long.high$type, y = EpidExt.long.high$Value, color = EpidExt.long.high$type)) + 
        scale_x_discrete(limits = levels(EpidExt.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        labs(x = "Run Type", y = "# Infected Counties", 
             title = eval(substitute(paste0("# Infected counties >",v), list(v=EpidExt_cutoff)))) + 
        scale_color_manual(values = cbPalette)
    })
    dev.off()
    
    EpidExt_Median_RepType=tapply(EpidExt.long$Value,EpidExt.long$type,FUN=median)
    EpidExt_Upper_RepType<-tapply(EpidExt.long$Value,EpidExt.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
    
    
    
    ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
    # two-column dataframes with fips for mapping
    
    EpidExt_medians=NULL
    EpidExt_uppers=NULL
    
    for (i in 1:length(run.types)){
      
      df <- EpidExt[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
      df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
      newname <- df1[1,4]
      df.new <- df[,-grep("type_", colnames(df))]
      
      # Calculate median and upper 2.5%
      df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
      df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
      
      # Create the vectors for map scales
      EpidExt_medians=c(EpidExt_medians,df.new$median)
      EpidExt_uppers=c(EpidExt_uppers,df.new$upper)
      
      # save datasets for mapping, one for each level with medians and uppers with name of control type
      assign(paste0("EpidExt_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
      assign(paste0("EpidExt_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
      
      rm(df, df1,df.new,newname)
    }
    
    ## Median maps
    
    # Create a vector for scale for median maps
    EpidExt_median_values <- round(c(0, 1, 5, mean(EpidExt_medians, na.rm = TRUE), 
                                     mean(EpidExt_medians, na.rm = TRUE) + sd(EpidExt_medians, na.rm = TRUE), 
                                     mean(EpidExt_medians, na.rm = TRUE) + 3*sd(EpidExt_medians, na.rm = TRUE), 
                                     mean(EpidExt_medians, na.rm = TRUE) + 6*sd(EpidExt_medians, na.rm = TRUE), 
                                     max(EpidExt_medians, na.rm = TRUE)), 0)
    EpidExt_median_values<-unique(sort(EpidExt_median_values))
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("EpidExt_med_",run.types[i])) 
      jpeg(paste0(path_output, paste0("EpidExt_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = EpidExt_median_values, color.sequence = palette, 
                  legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
      dev.off()
    }
    
    ## Upper 2.5% maps
    
    # Create a vector for scale for upper maps
    EpidExt_upper_values <- round(c(0, 1, 5, mean(EpidExt_uppers, na.rm = TRUE), 
                                    mean(EpidExt_uppers, na.rm = TRUE) + sd(EpidExt_uppers, na.rm = TRUE), 
                                    mean(EpidExt_uppers, na.rm = TRUE) + 3*sd(EpidExt_uppers, na.rm = TRUE), 
                                    mean(EpidExt_uppers, na.rm = TRUE) + 6*sd(EpidExt_uppers, na.rm = TRUE), # higher than the max
                                    max(EpidExt_uppers, na.rm = TRUE)), 0)
    if(6*sd(EpidExt_uppers, na.rm = TRUE)> max(EpidExt_uppers, na.rm = TRUE)) EpidExt_upper_values<-EpidExt_upper_values[c(1:6,8)]
    EpidExt_uppers<-unique(sort(EpidExt_upper_values))
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("EpidExt_upper_",run.types[i])) 
      jpeg(paste0(path_output, paste0("EpidExt_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = EpidExt_upper_values, color.sequence = palette, 
                  legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
      dev.off()
    }
    
    
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #####  Number Areas affected by movement ban (summary)      #######
  
  #####################################################################################################
  #####################################################################################################
  
  if (movementBan == TRUE){
    
    setwd(path_output)
    
    MB=county.summary[,grepl( "fips|polyname|shipBanImplemented|Type" , names( county.summary ) )]
    
    deleteme=c(0,rep(NA,ncol(MB)-1))
    
    # this gets rid of adjacent "type" columns, which indicate that a movement ban wasn't implemented
    for(i in 2:ncol(MB)){
      deleteme[i]=ifelse(grepl("Type", colnames(MB)[i-1]) & grepl("Type", colnames(MB)[i])|
                           grepl("Type", colnames(MB)[i]) & grepl("Type", colnames(MB)[i+1]),1,0)
    }
    
    MB=MB[,deleteme==0]
    
    if (ncol(MB)>2){
      
      names(MB)[seq(3,ncol(MB),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(MB))/2-1))))
      names(MB)[seq(4,ncol(MB),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(MB))/2-1))))
      if (export.datafiles == 1 | export.datafiles == 3) {write.csv(MB, "MovementBan_NumGeographyAffected.csv",row.names=F)}
    
      # Convert to long format for plotting
      MB.long <- reshape(MB, idvar = c("fips", "polyname"), direction = "long", 
                         v.names = c("type", "Value"), 
                         varying = list(c(grep("type_", colnames(MB))), 
                                        c(grep("run_", colnames(MB)))))
      
      rownames(MB.long) <- seq(1:nrow(MB.long))
      MB.long <- MB.long[,-3]
      MB.long <- MB.long[!is.na(MB.long$Value),] # lots of NAs
      if (export.datafiles == 2 | export.datafiles == 3) {write.csv(MB.long,"MovementBan_NumGeographyAffected_long.csv",row.names = F)}
      
      
      summary(MB.long$Value)
      table(MB.long$Value)
      hist(MB.long$Value)
      
      # Set Thresholds
        
      MB.long.overMin <- MB.long[MB.long$Value> MB_min,]
      MB.long.low <- MB.long[MB.long$Value <= MB_cutoff,]
      MB.long.high <- MB.long[MB.long$Value > MB_cutoff,]
      
      ## Histograms ##
      jpeg("MB_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(MB.long.high$Value, breaks = 25, xlab = "# Areas with Movement Ban", 
           main = paste0("# Movement Ban-affected geographies > ",as.character(MB_cutoff)))
      dev.off()
      
      jpeg("MB_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(MB.long.low$Value, breaks = 25, 
           xlab = "# Areas with Movement Ban", 
           main = paste0("# Movement Ban-affected geographies \u2264",as.character(MB_cutoff)))
      dev.off()
      
      jpeg("MB_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(MB.long.overMin$Value, breaks = 25, xlab = "# Areas with Movement Ban", 
           main = paste0("# Movement Bans >",as.character(MB_min)))
      dev.off()
      
      
      ## Violin plots ##
      jpeg("MB_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = MB.long.overMin, aes(x = MB.long.overMin$type, y = MB.long.overMin$Value, color = MB.long.overMin$type)) + 
          scale_x_discrete(limits = levels(MB.long.overMin$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Areas with Movement Ban", 
               title = eval(substitute(paste0("# Areas with Movement Ban >",v), list(v=MB_min)))) +  
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      jpeg("MB_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = MB.long.low, aes(x = MB.long.low$type, y = MB.long.low$Value, color = MB.long.low$type)) + 
          scale_x_discrete(limits = levels(MB.long.low$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Areas with Movement Ban", title = paste0("# Areas with Movement Ban \u2264", as.character(MB_cutoff) )) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      jpeg("MB_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = MB.long.high, aes(x = MB.long.high$type, y = MB.long.high$Value, color = MB.long.high$type)) + 
          scale_x_discrete(limits = levels(MB.long.high$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Areas with Movement Ban", 
               title = eval(substitute(paste0("# Areas with Movement Bans >",v), list(v=MB_cutoff)))) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      
      
      # Calculate the median of the high and low datasets for each control type for 
      # summary table (not used elsewhere)
      MB.run.types=run.types[grep("MvmtBan",run.types)] 
      
      
      MB_Median_RepType=tapply(MB.long$Value,MB.long$type,FUN=median)
      MB_Upper_RepType<-tapply(MB.long$Value,MB.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
      
      
      ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
      # two-column dataframes with fips for mapping
      
      MB_medians=NULL
      MB_uppers=NULL
      
      for (i in 1:length(MB.run.types)){
        
        df <- MB[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
        df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
        newname <- df1[1,4]
        df.new <- df[,-grep("type_", colnames(df))]
        
        # Calculate median and upper 2.5%
        df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
        df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
        
        # Create the vectors for map scales
        MB_medians=c(MB_medians,df.new$median)
        MB_uppers=c(MB_uppers,df.new$upper)
        
        # save datasets for mapping, one for each level with medians and uppers with name of control type
        assign(paste0("MB_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
        assign(paste0("MB_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
        
        rm(df, df1,df.new,newname)
      }
      
      ## Median maps
      
      # Create a vector for scale for median maps
      MB_median_values <- round(c(0, 1, 10, mean(MB_medians, na.rm = TRUE), 
                                  mean(MB_medians, na.rm = TRUE) + sd(MB_medians, na.rm = TRUE), 
                                  mean(MB_medians, na.rm = TRUE) + 3*sd(MB_medians, na.rm = TRUE), 
                                  mean(MB_medians, na.rm = TRUE) + 6*sd(MB_medians, na.rm = TRUE), 
                                  max(MB_medians, na.rm = TRUE)), 0)
      MB_median_values<-sort(unique(MB_median_values))
      
      for (i in 1:length(MB.run.types)){
        name_df=get(paste0("MB_med_",MB.run.types[i])) 
        if (var(name_df$median, na.rm=TRUE)==0) {name_df$median[1]<-(name_df$median[1]+0.001)} # map_by_fips doesn't work with all zeros, so add a tiny bit. 
        jpeg(paste0(path_output, paste0("MB_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = MB_median_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
        dev.off()
      }
      
      ## Upper 2.5% maps
      
      # Create a vector for scale for median maps
      MB_upper_values <- round(c(0, 1, 10, mean(MB_uppers, na.rm = TRUE), 
                                 mean(MB_uppers, na.rm = TRUE) + sd(MB_uppers, na.rm = TRUE), 
                                 mean(MB_uppers, na.rm = TRUE) + 3*sd(MB_uppers, na.rm = TRUE), 
                                 mean(MB_uppers, na.rm = TRUE) + 6*sd(MB_uppers, na.rm = TRUE), 
                                 max(MB_uppers, na.rm = TRUE)), 0)
      if(6*sd(MB_uppers, na.rm = TRUE)> max(MB_uppers, na.rm = TRUE)) MB_upper_values<-MB_upper_values[c(1:6,8)]
      MB_upper_values<-sort(MB_upper_values)
      
      for (i in 1:length(MB.run.types)){
        name_df=get(paste0("MB_upper_",MB.run.types[i])) 
        jpeg(paste0(path_output, paste0("MB_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = MB_upper_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
        dev.off()
      }
    }
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #####  Premises Culled (summary)      #######
  
  #####################################################################################################
  #####################################################################################################
  
  # Problems (fixed, see comments)
  if (premisesCulled == TRUE) {
    setwd(path_output)
    # Filtering out runs without culling first
    
    PremCull=county.summary[,grepl( "fips|polyname|cullImplemented|Type" , names( county.summary ) )]
    
    PremCull=PremCull[,!grepl( "DCSubset" , names( PremCull ) )]
    
    
    deleteme=c(0,rep(NA,ncol(PremCull)-1))
    
    # this gets rid of adjacent "type" columns, which indicate that a cull wasn't implemented
    for(i in 2:ncol(PremCull)){
      deleteme[i]=ifelse(grepl("Type", colnames(PremCull)[i-1]) & grepl("Type", colnames(PremCull)[i])|
                           grepl("Type", colnames(PremCull)[i]) & grepl("Type", colnames(PremCull)[i+1]),1,0)
    }
    
    PremCull=PremCull[,deleteme==0]
    
    if(ncol(PremCull)>2){
      
      names(PremCull)[seq(3,ncol(PremCull),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(PremCull))/2-1))))
      names(PremCull)[seq(4,ncol(PremCull),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(PremCull))/2-1))))
      if (export.datafiles == 1 | export.datafiles == 3) {write.csv(PremCull,"PremisesCulled.csv",row.names=F)}

      
      PremCull.long <- reshape(PremCull, idvar = c("fips", "polyname"), direction = "long", 
                               v.names = c("type", "Value"), 
                               varying = list(c(grep("type_", colnames(PremCull))), 
                                              c(grep("run_", colnames(PremCull)))))
      
      rownames(PremCull.long) <- seq(1:nrow(PremCull.long))
      PremCull.long <- PremCull.long[,-3]
      PremCull.long <- PremCull.long[!is.na(PremCull.long$Value),] 
      if (export.datafiles == 2 | export.datafiles == 3) {write.csv(PremCull.long,"PremisesCulled_long.csv",row.names = F)}

      
      summary(PremCull.long$Value)
      hist(PremCull.long$Value)
      
      ## Divide the data into subsets for plotting ## 
      
      # Set Thresholds
        
      PremCull.long.overMin <- PremCull.long[PremCull.long$Value > PremCull_min,]
      PremCull.long.low <- PremCull.long[PremCull.long$Value <= PremCull_cutoff,]
      PremCull.long.high <- PremCull.long[PremCull.long$Value > PremCull_cutoff,]
      
      # Histograms
      jpeg("PremCull_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(PremCull.long.high$Value, breaks = 25, xlab = "# Premises culled", 
           main = paste0("# Premises culled > ",as.character(PremCull_cutoff)))
      dev.off()
      
      jpeg("PremCull_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(PremCull.long.low$Value, breaks = 25, 
           xlab = "# Premises culled", 
           main = paste0("# Premises culled \u2264 ",as.character(PremCull_cutoff)))
      dev.off()
      
      jpeg("PremCull_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(PremCull.long.overMin$Value, breaks = 25, xlab = "# Premises culled", 
           main = paste0("# Premises culled > ",as.character(PremCull_min)))
      dev.off()
      
      # violin plot testing
      jpeg("PremCull_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = PremCull.long.overMin, aes(x = PremCull.long.overMin$type, y = PremCull.long.overMin$Value, color = PremCull.long.overMin$type)) + 
          scale_x_discrete(limits = levels(PremCull.long.overMin$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Culled Premises", 
               title = eval(substitute(paste0("# Culled premises >",v), list(v=PremCull_min)))) +  
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      jpeg("PremCull_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = PremCull.long.low, aes(x = PremCull.long.low$type, y = PremCull.long.low$Value, color = PremCull.long.low$type)) + 
          scale_x_discrete(limits = levels(PremCull.long.low$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Culled Premises", title = paste0("# Culled premises \u2264", as.character(PremCull_cutoff) )) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      jpeg("PremCull_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = PremCull.long.high, aes(x = PremCull.long.high$type, y = PremCull.long.high$Value, color = PremCull.long.high$type)) + 
          scale_x_discrete(limits = levels(PremCull.long.high$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Culled Premises", 
               title = eval(substitute(paste0("# Culled premises >",v), list(v=PremCull_cutoff)))) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      
      
      # Calculate the median of the high and low datasets for each control type for 
      # summary table (not used elsewhere)
    
      PremCull.run.types=run.types[grep("cull",run.types)]  
      
      
      PremCull_Median_RepType=tapply(PremCull.long$Value,PremCull.long$type,FUN=median)
      PremCull_Upper_RepType<-tapply(PremCull.long$Value,PremCull.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
      
      
      
      ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
      # two-column dataframes with fips for mapping
      
      PremCull_medians=NULL
      PremCull_uppers=NULL
      
      for (i in 1:length(PremCull.run.types)){
        df <- PremCull[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
        df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
        newname <- df1[1,4]
        df.new <- df[,-grep("type_", colnames(df))]
        
        # Calculate median and upper 2.5%
        df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
        df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
        
        # Create the vectors for map scales
        PremCull_medians=c(PremCull_medians,df.new$median)
        PremCull_uppers=c(PremCull_uppers,df.new$upper)
        
        # save datasets for mapping, one for each level with medians and uppers with name of control type
        assign(paste0("PremCull_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
        assign(paste0("PremCull_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
        
        rm(df, df1,df.new,newname)
      }
      
      ## Median maps
      
      # Create a vector for scale for median maps
      PremCull_median_values <- round(c(0, 1, 10, mean(PremCull_medians, na.rm = TRUE), 
                                        mean(PremCull_medians, na.rm = TRUE) + sd(PremCull_medians, na.rm = TRUE), 
                                        mean(PremCull_medians, na.rm = TRUE) + 3*sd(PremCull_medians, na.rm = TRUE), 
                                        mean(PremCull_medians, na.rm = TRUE) + 6*sd(PremCull_medians, na.rm = TRUE), 
                                        max(PremCull_medians, na.rm = TRUE)), 0)
      PremCull_median_values<-sort(unique(PremCull_median_values))
      
      for (i in 1:length(PremCull.run.types)){
        name_df=get(paste0("PremCull_med_",PremCull.run.types[i]))
        if (var(name_df$median, na.rm=TRUE)==0) {name_df$median[1]<-(name_df$median[1]+0.001)} # map_by_fips doesn't work with all zeros, so add a tiny bit.
        jpeg(paste0(path_output, paste0("PremCull_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = PremCull_median_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
        dev.off()
      }
      
      ## Upper 2.5% maps
      
      # Create a vector for scale for upper maps
      PremCull_upper_values <- round(c(0, 1, 10, mean(PremCull_uppers, na.rm = TRUE), 
                                       mean(PremCull_uppers, na.rm = TRUE) + sd(PremCull_uppers, na.rm = TRUE), 
                                       mean(PremCull_uppers, na.rm = TRUE) + 3*sd(PremCull_uppers, na.rm = TRUE), 
                                       mean(PremCull_uppers, na.rm = TRUE) + 6*sd(PremCull_uppers, na.rm = TRUE), # higher than the max
                                       max(PremCull_uppers, na.rm = TRUE)), 0)
      if(6*sd(PremCull_uppers, na.rm = TRUE)> max(PremCull_uppers, na.rm = TRUE)) PremCull_upper_values<-PremCull_upper_values[c(1:6,8)]
      PremCull_upper_values<-sort(PremCull_upper_values)
      
      for (i in 1:length(PremCull.run.types)){
        name_df=get(paste0("PremCull_upper_",PremCull.run.types[i])) 
        jpeg(paste0(path_output, paste0("PremCull_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = PremCull_upper_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
        dev.off()
      }
    } 
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #####  Premises Vaccinated (summary)      #######
  
  #####################################################################################################
  #####################################################################################################
  
  if (premisesVax == TRUE){
    # Filtering out runs without vaccination first
    
    PremVax=county.summary[,grepl( "fips|polyname|vaxImplemented|Type" , names( county.summary ) )]
    PremVax=PremVax[,!grepl( "DCSubset" , names( PremVax ) )]
 
    
    deleteme=c(0,rep(NA,ncol(PremVax)-1))
    
    # this gets rid of adjacent "type" columns, which indicate that a cull wasn't implemented
    for(i in 2:ncol(PremVax)){
      deleteme[i]=ifelse(grepl("Type", colnames(PremVax)[i-1]) & grepl("Type", colnames(PremVax)[i])|
                           grepl("Type", colnames(PremVax)[i]) & grepl("Type", colnames(PremVax)[i+1]),1,0)
    }
    
    PremVax=PremVax[,deleteme==0]
    
    setwd(path_output)
    
    
    if(ncol(PremVax)>2){
      
      names(PremVax)[seq(3,ncol(PremVax),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(PremVax))/2-1))))
      names(PremVax)[seq(4,ncol(PremVax),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(PremVax))/2-1))))
      if (export.datafiles == 1 | export.datafiles == 3) {write.csv(PremVax,"PremisesVaccinated.csv",row.names=F)}

      PremVax.long <- reshape(PremVax, idvar = c("fips", "polyname"), direction = "long", 
                              v.names = c("type", "Value"), 
                              varying = list(c(grep("type_", colnames(PremVax))), 
                                             c(grep("run_", colnames(PremVax)))))
      
      rownames(PremVax.long) <- seq(1:nrow(PremVax.long))
      PremVax.long <- PremVax.long[,-3]
      PremVax.long <- PremVax.long[!is.na(PremVax.long$Value),] 
      if (export.datafiles == 2 | export.datafiles == 3) {write.csv(PremVax.long,"PremisesVaccinated_long.csv",row.names = F)}

      summary(PremVax.long$Value)
      hist(PremVax.long$Value)
      
      ## Divide the data into subsets for plotting ## 
      
      PremVax.long.overMin <- PremVax.long[PremVax.long$Value> PremVax_min,]
      PremVax.long.low <- PremVax.long[PremVax.long$Value <= PremVax_cutoff,]
      PremVax.long.high <- PremVax.long[PremVax.long$Value > PremVax_cutoff,]
      
      ## Histograms ##
      jpeg("PremVax_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(PremVax.long.high$Value, breaks = 25, xlab = "# Premises vaccinated", 
           main = paste0("# Premises vaccinated > ",as.character(PremVax_cutoff)))
      dev.off()
      
      jpeg("PremVax_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(PremVax.long.low$Value, breaks = 25, 
           xlab = "# Premises vaccinated", 
           main = paste0("# Premises vaccinated \u2264 ",as.character(PremVax_cutoff)))
      dev.off()
      
      jpeg("PremVax_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
      hist(PremVax.long$Value[PremVax.long$Value > PremVax_min], breaks = 25, xlab = "# Premises vaccinated", 
           main = paste0("# Premises vaccinated > ",as.character(PremVax_min)))
      dev.off()
      
      
      ## Violin plots ##
      jpeg("PremVax_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = PremVax.long.overMin, aes(x = PremVax.long.overMin$type, y = PremVax.long.overMin$Value, color = PremVax.long.overMin$type)) + 
          scale_x_discrete(limits = levels(PremVax.long.overMin$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Vaccinated Premises", 
               title = eval(substitute(paste0("# Vaccinated premises >",v), list(v=PremVax_min)))) +  
          scale_color_manual(values = cbPalette) 
      })
      dev.off()
      
      jpeg("PremVax_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = PremVax.long.low, aes(x = PremVax.long.low$type, y = PremVax.long.low$Value, color = PremVax.long.low$type)) + 
          scale_x_discrete(limits = levels(PremVax.long.low$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Vaccinated Premises", 
               title = paste0("# Vaccinated premises \u2264", as.character(PremVax_cutoff) )) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      jpeg("PremVax_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
      print({
        ggplot() + geom_violin(data = PremVax.long.high, aes(x = PremVax.long.high$type, y = PremVax.long.high$Value, color = PremVax.long.high$type)) + 
          scale_x_discrete(limits = levels(PremVax.long.high$type)) + theme_bw() + 
          theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
          labs(x = "Run Type", y = "# Vaccinated Premises", 
               title = eval(substitute(paste0("# Vaccinated premises >",v), list(v=PremVax_cutoff)))) + 
          scale_color_manual(values = cbPalette)
      })
      dev.off()
      
      
      
      # Calculate the median of the high and low datasets for each control type for 
      # summary table (not used elsewhere)
      PremVax.run.types=run.types[grep("VAX",run.types)] 
      PremVax.run.types=run.types[grep("Vax",run.types)] 
      
      
      PremVax_Median_RepType=tapply(PremVax.long$Value,PremVax.long$type,FUN=median)
      PremVax_Upper_RepType<-tapply(PremVax.long$Value,PremVax.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
      
      
      ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
      # two-column dataframes with fips for mapping
      
      PremVax_medians=NULL
      PremVax_uppers=NULL
      
      for (i in 1:length(PremVax.run.types)){
        
        df <- PremVax[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
        df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
        newname <- df1[1,4]
        df.new <- df[,-grep("type_", colnames(df))]
        
        # Calculate median and upper 2.5%
        df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
        df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
        
        # Create the vectors for map scales
        PremVax_medians=c(PremVax_medians,df.new$median)
        PremVax_uppers=c(PremVax_uppers,df.new$upper)
        
        # save datasets for mapping, one for each level with medians and uppers with name of control type
        assign(paste0( "PremVax_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
        assign(paste0( "PremVax_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
        
        rm(df, df1,df.new,newname)
      }
      
      ## Median maps
      
      # Create a vector for scale for median maps
      PremVax_median_values <- round(c(0, 1, 10, mean(PremVax_medians, na.rm = TRUE), 
                                       mean(PremVax_medians, na.rm = TRUE) + sd(PremVax_medians, na.rm = TRUE), 
                                       mean(PremVax_medians, na.rm = TRUE) + 3*sd(PremVax_medians, na.rm = TRUE), 
                                       mean(PremVax_medians, na.rm = TRUE) + 6*sd(PremVax_medians, na.rm = TRUE), 
                                       max(PremVax_medians, na.rm = TRUE)), 0)
      PremVax_median_values<-sort(unique(PremVax_median_values))
      
      for (i in 1:length(PremVax.run.types)){
        name_df=get(paste0("PremVax_med_",PremVax.run.types[i]))
        if (var(name_df$median, na.rm=TRUE)==0) {name_df$median[1]<-(name_df$median[1]+0.001)} # map_by_fips doesn't work with all zeros, so add a tiny bit.
        jpeg(paste0(path_output, paste0("PremVax_Median_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = PremVax_median_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.5, legend.width = 1)
        dev.off()
      }
      
      ## Upper 2.5% maps
      
      # Create a vector for scale for upper maps
      PremVax_upper_values <- round(c(0, 1, 10, mean(PremVax_uppers, na.rm = TRUE), 
                                      mean(PremVax_uppers, na.rm = TRUE) + sd(PremVax_uppers, na.rm = TRUE), 
                                      mean(PremVax_uppers, na.rm = TRUE) + 3*sd(PremVax_uppers, na.rm = TRUE), 
                                      mean(PremVax_uppers, na.rm = TRUE) + 6*sd(PremVax_uppers, na.rm = TRUE), 
                                      max(PremVax_uppers, na.rm = TRUE)), 0)
      if(6*sd(PremVax_uppers, na.rm = TRUE)> max(PremVax_uppers, na.rm = TRUE)) PremVax_upper_values<-PremVax_upper_values[c(1:6,8)]
      PremVax_upper_values<-sort(PremVax_upper_values)
      
      for (i in 1:length(PremVax.run.types)){
        name_df=get(paste0("PremVax_upper_",PremVax.run.types[i])) 
        jpeg(paste0(path_output, paste0("PremVax_Upper_Map_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
        setwd(path0)
        map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                    missing.include = TRUE, color.break.type = "values", 
                    color.break.values = PremVax_upper_values, color.sequence = palette, 
                    legend.spacing = 4.5, legend.shrink = 0.3, legend.width = 1)
        dev.off()
      }
    }
  }

  
  #####################################################################################################
  #####################################################################################################
  
  #####  List detail files and load FLAPS      #######
  
  #####################################################################################################
  #####################################################################################################
  
  if(animalsInfected == TRUE | countyRisk == TRUE | localSpread == TRUE){
    
    ## Find all the detail files in the Files_To_Process" directory ##
    detail.fnames <- list.files(path = pathfiles, recursive = TRUE, pattern = "_detail.txt", full.names = FALSE)
    
    ## Read in and format all FLAPS files ##
    # These can be merged with detail files to identify characteristics of (potentially) exposed/infected farms. 
    
    setwd(path0)
    
    for(flap in 1:10){
      flname=paste0("f",flap)
      assign(flname,fread(paste0("FLAPS/FLAPS12_Quarterly_USDOS_format_",sprintf("%04d",flap),".txt"))) 
    }
    
    for (f in 1:10){
      df=get(paste0("f",f))
      df$anim=df$b_Q3+df$d_Q3
      df$County_fips[which(df$County_fips==46113)]<- 46102
      assign(paste0("f",f),df)
    }
    rm(df)
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #####  Number of Animals Infected (Combines summary and Detail Files)      #######
  # End result is one line per rep, so treat like a summary metric. 
  
  #####################################################################################################
  #####################################################################################################
  
  if(animalsInfected == TRUE){
    
    # Initiate the county.anim file with the first single file 
    # read the summary file
    setwd(pathfiles)
    sum.file <- fread(summary.files[1], header = TRUE, select = c("Rep","Num_Inf", "Seed_Farms", "Seed_FIPS"))
    
    setwd(path0)
    
    # Select the right FLAPS file
    FLAP <- substr(unlist(strsplit(summary.files[1],"format_"))[2],1,4)  
    
    flaps <- get(paste0("f",as.numeric(FLAP)))
    
    # merge FLAPS with summary file to get # animals on seed farm
    sum.flaps <- merge(sum.file,flaps[,c("Id", "anim")],by.x="Seed_Farms",by.y="Id",all.x=TRUE)
    
    # Read detail file
    setwd(pathfiles)
    det.file <- fread(detail.fnames[1], header = TRUE, select = c("Rep","ExposedID", "ExposedCounty", "SourceID", "SourceCounty", "ControlPrevented"))
    
    # Filter on not prevented and unique ExposureID by Rep
    det.file <- det.file[det.file$ControlPrevented == "none",]
    det.file <- unique(as.data.table(det.file),by=c("Rep","ExposedID")) 
    
    # Summarize the number of animals infected per Rep
    rep.totals <- merge(flaps[,c("Id", "County_fips", "anim")], det.file[,c("Rep","ExposedID", "ExposedCounty")], 
                        by.x = c("Id", "County_fips"), by.y = c("ExposedID", "ExposedCounty"), all = TRUE)
    rep.totals <- as.data.frame(aggregate(rep.totals$anim, by = list(rep.totals$Rep), FUN = sum))
    colnames(rep.totals) <- c("Rep", "anim")
    
    # If a Rep's Num_Inf >1, get the additional (non seed farm) animals infected in that rep from the detail file
    for(row in 1:max(sum.flaps$Rep)){
      if(sum.flaps$Num_Inf[row]>1){sum.flaps$anim[row]=sum.flaps$anim[row]+rep.totals$anim[which(rep.totals$Rep==sum.flaps$Rep[row])]}
    }
    
    Anim=sum.flaps[,c("Seed_FIPS","anim")]
    
    # Get type
    Anim$Type <- unlist(strsplit(summary.files[1], "_FLAPS"))[1]
    Anim<-as.data.frame(Anim)
    
    
    # Read through each summary/detail file 
    for (file in 2:length(detail.fnames)){
      # read summary file
      sum.file <- fread(summary.files[file], header = TRUE, select = c("Rep","Num_Inf", "Seed_Farms", "Seed_FIPS"))
      
      # Select the right FLAPS file
      FLAP <- substr(unlist(strsplit(summary.files[file],"format_"))[2],1,4)  
      flaps = get(paste0("f",as.numeric(FLAP)))
      
      # merge FLAPS with summary file 
      sum.flaps=merge(sum.file,flaps[,c("Id", "anim")],by.x="Seed_Farms",by.y="Id",all.x=TRUE)
      
      # Read detail file
      det.file <- fread(detail.fnames[file], header = TRUE, select = c("Rep","ExposedID", "ExposedCounty", "SourceID", "SourceCounty", "ControlPrevented"))
      
      # Filter on not prevented and unique ExposureID by Rep
      det.file <- det.file[det.file$ControlPrevented == "none",]
      det.file <- unique(as.data.table(det.file),by=c("Rep","ExposedID")) 
      
      # Summarize the number of animals infected per Rep
      rep.totals <- merge(flaps[,c("Id", "County_fips", "anim")], det.file[,c("Rep","ExposedID", "ExposedCounty")], 
                          by.x = c("Id", "County_fips"), by.y = c("ExposedID", "ExposedCounty"), all = TRUE)
      rep.totals <- as.data.frame(aggregate(rep.totals$anim, by = list(rep.totals$Rep), FUN = sum))
      colnames(rep.totals) <- c("Rep", "anim")
      
      # If a Rep's Num_Inf >1, get the additional (non seed farm) animals infected in that rep from the detail file
      for(row in 1:max(sum.flaps$Rep)){
        if(sum.flaps$Num_Inf[row]>1){sum.flaps$anim[row]=sum.flaps$anim[row]+rep.totals$anim[which(rep.totals$Rep==sum.flaps$Rep[row])]}
      }
      
      run.anim=sum.flaps[,c("Seed_FIPS","anim")]
      
      # Get type
      run.anim$Type <- unlist(strsplit(summary.files[file], "_FLAPS"))[1]
      run.anim <- as.data.frame(run.anim)
      
      Anim=merge(Anim,run.anim, by="Seed_FIPS")
      
    } # End summary/detail file loop
    
    
    # now it's generally shaped like the summary file metrics. Needs the below, and is also 3049 rows versus 3086
    # rename column 1
    names(Anim)[1]="fips"
    # Add col 2 of polynames 
    
    Anim=merge(county.summary[,1:2],Anim,by="fips",all=TRUE) # This brings the total to the same as the other metrics. This and other metrics have NA values for 27 counties with no premises 
    Anim <- unique(as.data.table(Anim),by=c("fips"))  # cut to a single entry for each FIPS 
    Anim=as.data.frame(Anim)
    
    names(Anim)[seq(3,ncol(Anim),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(Anim))/2-1))))
    names(Anim)[seq(4,ncol(Anim),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(Anim))/2-1))))
    
    setwd(path_output)
    
    
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(Anim, "NumAnimalsInfected.csv",row.names=F)}

    Anim.long <- reshape(Anim, idvar = c("fips", "polyname"), direction = "long", 
                         v.names = c("type", "Value"), 
                         varying = list(c(grep("type_", colnames(Anim))), 
                                        c(grep("run_", colnames(Anim)))))
    
    rownames(Anim.long) <- seq(1:nrow(Anim.long))
    Anim.long <- Anim.long[,-3]
    Anim.long <- Anim.long[!is.na(Anim.long$Value),] 
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(Anim.long,"AnimalsInfected_long.csv",row.names = F)}


    Anim.long.overMin <- Anim.long[Anim.long$Value> Anim_min,]
    Anim.long.low <- Anim.long[Anim.long$Value <= Anim_cutoff,]
    Anim.long.high <- Anim.long[Anim.long$Value > Anim_cutoff,]
    
    # Histograms
    jpeg("Anim_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(Anim.long.high$Value, breaks = 25, xlab = "# Animals Infected", 
         main = paste0("# Animals Infected > ",as.character(Anim_cutoff)))
    dev.off()
    
    jpeg("Anim_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(Anim.long.low$Value, breaks = 25, 
         xlab = "# Animals infected", 
         main = paste0("# Animals infected \u2264",as.character(Anim_cutoff)))
    
    dev.off()
    # 
    jpeg("Anim_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(Anim.long$Value[Anim.long$Value > 0], breaks = 25, xlab = "# Animals infected",
         main = paste0("# Animals infected > ",Anim_min))
    dev.off()
    
    
    # violin plots
    jpeg("Anim_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = Anim.long.overMin, aes(x = Anim.long.overMin$type, y = Anim.long.overMin$Value, 
                                                           color = Anim.long.overMin$type)) +
        scale_x_discrete(limits = levels(Anim.long.overMin$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        labs(x = NULL, y = "# Animals infected", title = paste0("# Animals infected > ",Anim_min)) + scale_color_manual(values = cbPalette)
    })
    dev.off()
    
    jpeg("Anim_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = Anim.long.low, aes(x = Anim.long.low$type, y = Anim.long.low$Value, color = Anim.long.low$type)) + 
        scale_x_discrete(limits = levels(Anim.long.low$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        scale_color_manual(values = cbPalette) + labs(x = NULL, y = "# Infected Animals", 
                                                      title = paste0("# Animals infected  \u2264 ",as.character(Anim_cutoff)))
    })
    dev.off()
    
    jpeg("Anim_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = Anim.long.high, aes(x = Anim.long.high$type, y = Anim.long.high$Value, color = Anim.long.high$type)) + 
        scale_x_discrete(limits = levels(Anim.long.high$type)) + theme_bw() + 
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
        scale_color_manual(values = cbPalette) + labs(x = NULL, y = "# Infected Animals", 
                                                      title = paste0("# Animals infected > ",as.character(Anim_cutoff)))
    })
    dev.off()
    
    

    Anim_Median_RepType=tapply(Anim.long$Value,Anim.long$type,FUN=median)
    Anim_Upper_RepType<-tapply(Anim.long$Value,Anim.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
    
    
    
    ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into 
    # two-column dataframes with fips for mapping
    
    Anim_medians=NULL
    Anim_uppers=NULL
    
    for (i in 1:length(run.types)){
      df <- Anim[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
      df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
      newname <- df1[1,4]
      df.new <- df[,-grep("type_", colnames(df))]
      
      # Calculate median and upper 2.5%
      df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
      df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
      
      # Create the vectors for map scales
      Anim_medians=c(Anim_medians,df.new$median)
      Anim_uppers=c(Anim_uppers,df.new$upper)
      
      # save datasets for mapping, one for each level with medians and uppers with name of control type
      assign(paste0("Anim_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
      assign(paste0("Anim_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
      
      rm(df, df1,df.new,newname)
    }
    
    ## Median maps
    
    # Create a vector for scale for median maps
    Anim_median_values <- round(c(0, 100, 1000, mean(Anim_medians, na.rm = TRUE), 
                                  mean(Anim_medians, na.rm = TRUE) + sd(Anim_medians, na.rm = TRUE), 
                                  mean(Anim_medians, na.rm = TRUE) + 3*sd(Anim_medians, na.rm = TRUE), 
                                  mean(Anim_medians, na.rm = TRUE) + 6*sd(Anim_medians, na.rm = TRUE), 
                                  max(Anim_medians, na.rm = TRUE)), 0)
    if(6*sd(Anim_median_values, na.rm = TRUE)> max(Anim_median_values, na.rm = TRUE)) {Anim_median_values<-Anim_median_values[c(1:6,8)]}
    if(3*sd(Anim_median_values, na.rm = TRUE)> max(Anim_median_values, na.rm = TRUE)) {Anim_median_values<-Anim_median_values[c(1:5,7)]}
    Anim_median_values<-sort(unique(Anim_median_values))
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("Anim_med_",run.types[i])) 
      jpeg(paste0(path_output, paste0("Anim_Median_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = Anim_median_values, color.sequence = palette, 
                  legend.spacing = 6.5, legend.shrink = 0.4, legend.width = 1, cex.axis = 0.8)
      dev.off()
    }
    
    ## Upper 2.5% maps
    
    # Create a vector for scale for median maps
    Anim_upper_values <- round(c(0, 100, 1000, 100000,mean(Anim_uppers, na.rm = TRUE), 
                                 mean(Anim_uppers, na.rm = TRUE) + sd(Anim_uppers, na.rm = TRUE), 
                                 mean(Anim_uppers, na.rm = TRUE) + 3*sd(Anim_uppers, na.rm = TRUE), 
                                 mean(Anim_uppers, na.rm = TRUE) + 6*sd(Anim_uppers, na.rm = TRUE), 
                                 max(Anim_uppers, na.rm = TRUE)+1), 0)
    if(6*sd(Anim_upper_values, na.rm = TRUE)> max(Anim_upper_values, na.rm = TRUE)) {Anim_upper_values<-Anim_upper_values[c(1:7,9)]}
    if(3*sd(Anim_upper_values, na.rm = TRUE)> max(Anim_upper_values, na.rm = TRUE)) {Anim_upper_values<-Anim_upper_values[c(1:6,8)]}
    Anim_upper_values<-sort(Anim_upper_values)
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("Anim_upper_",run.types[i])) 
      jpeg(paste0(path_output, paste0("Anim_Up_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30", 
                  missing.include = TRUE, color.break.type = "values", 
                  color.break.values = Anim_upper_values, color.sequence = palette, 
                  legend.spacing = 6.5, legend.shrink = 0.4, legend.width = 1, cex.axis = 0.8)
      dev.off()
    }
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #####  County Risk (Detail & Summary files)      #######
  
  #####################################################################################################
  #####################################################################################################
  
  if (countyRisk == TRUE){
    setwd(pathfiles)
    
    # read the summary file
    sum.file <- as.data.frame(fread(summary.files[1], header = TRUE, select = c("Rep", "Seed_FIPS")))
    
    # read the detail file 
    det.file <- as.data.frame(fread(detail.fnames[1], header = TRUE, select = c("Rep","SourceCounty","ExposedCounty", "ControlPrevented")))
    
    # remove prevented exposures from the detail file 
    det.file <- det.file[det.file$ControlPrevented == "none",] # successful exposures
    
    # Merge in the seed FIPS from the summary file
    det.file <- merge(sum.file,det.file,by="Rep")
    
    # exclude any with exposed county = seed county
    det.file <- det.file[det.file$Seed_FIPS!=det.file$ExposedCounty,]
    
    # Get unique exposures per rep (only count an exposed county once per rep)
    CountyRisk <- unique(as.data.table(det.file),by=c("Rep","ExposedCounty")) 
    
    # table it to get the total number of seeded reps in which a county was exposed (exlucding reps in which that county was the seed)
    CountyRisk <- as.data.frame(table(CountyRisk$ExposedCounty))
    
    # merge with county list to retrieve counties that are missing
    CountyRisk <- merge(sum.file, CountyRisk, by.x = "Seed_FIPS", by.y="Var1",all=TRUE )
    
    # Replace Frequency valeus of NA (i.e. counties that were not exposed in any reps) with 0
    CountyRisk$Freq[is.na(CountyRisk$Freq)] <- 0
    
    # calculate county risk as a proportion by dividing the number of reps in which it was exposed by the number of reps -1 (for the rep in which it was the seed)
    CountyRisk$CountyRisk <- CountyRisk$Freq/(length(CountyRisk$Seed_FIPS)-1)
    
    # format 
    CountyRisk <- CountyRisk[,c("Seed_FIPS","CountyRisk")]
    CountyRisk$Type <- unlist(strsplit(summary.files[1], "_FLAPS"))[1]
    
    
    for(i in 2:length(detail.fnames)) {
      sum.file <- as.data.frame(fread(summary.files[i], header = TRUE, select = c("Rep", "Seed_FIPS")))
      det.file <- as.data.frame(fread(detail.fnames[i], header = TRUE, select = c("Rep","SourceCounty","ExposedCounty", "ControlPrevented")))
      det.file <- det.file[det.file$ControlPrevented == "none",] # successful exposures
      det.file <- merge(sum.file,det.file,by="Rep")
      det.file <- det.file[det.file$Seed_FIPS!=det.file$ExposedCounty,]
      cr <- unique(as.data.table(det.file),by=c("Rep","ExposedCounty")) 
      cr <- as.data.frame(table(cr$ExposedCounty))
      cr <- merge(sum.file, cr, by.x = "Seed_FIPS", by.y="Var1",all=TRUE )
      cr$Freq[is.na(cr$Freq)] <- 0
      cr$CountyRisk <- cr$Freq/(length(cr$Seed_FIPS)-1) 
      cr <- cr[,c("Seed_FIPS","CountyRisk")]
      cr$Type <- unlist(strsplit(summary.files[i], "_FLAPS"))[1]
      
      # merge with previous iterations
      CountyRisk <- merge(CountyRisk,cr,by="Seed_FIPS",all=TRUE)
      
    }
    
    
    # Format like summary file metrics
    
    names(CountyRisk)[1]="fips"
    
    CountyRisk=merge(county.summary[,1:2],CountyRisk,by="fips",all=TRUE)
  
    names(CountyRisk)[seq(4,ncol(CountyRisk),2)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(CountyRisk))/2-1))))
    names(CountyRisk)[seq(3,ncol(CountyRisk),2)] <- paste0("run_", sprintf("%03d", seq(1:((ncol(CountyRisk))/2-1))))
    
    setwd(path_output)
    
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(CountyRisk, "CountyRisk.csv",row.names = F)}
    
    CountyRisk.long <- reshape(CountyRisk, idvar = c("fips", "polyname"), direction = "long",
                               v.names = c("type", "Value"),
                               varying = list(c(grep("type_", colnames(CountyRisk))),
                                              c(grep("run_", colnames(CountyRisk)))))
    
    rownames(CountyRisk.long) <- seq(1:nrow(CountyRisk.long))
    CountyRisk.long <- CountyRisk.long[,-3]
    CountyRisk.long <- CountyRisk.long[!is.na(CountyRisk.long$Value),] 
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(CountyRisk.long,"CountyRisk_long.csv",row.names = F)}


    CountyRisk.long.overMin <- CountyRisk.long[CountyRisk.long$Value > CountyRisk_min,]
    CountyRisk.long.low <- CountyRisk.long[CountyRisk.long$Value <= CountyRisk_cutoff,]
    CountyRisk.long.high <- CountyRisk.long[CountyRisk.long$Value > CountyRisk_cutoff,]
    
    # Histograms
    jpeg("CountyRisk_High_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(CountyRisk.long.high$Value, breaks = 25, xlab = "County Risk",
         main = paste0("County Risk > ",as.character(CountyRisk_cutoff)))
    dev.off()
    
    jpeg("CountyRisk_Low_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(CountyRisk.long.low$Value, breaks = min(25,length(unique(CountyRisk.long.low$Value))),
         xlab = "CountyRisk", main = paste0("County Risk \u2264 ",as.character(CountyRisk_cutoff)))
    dev.off()
    
    jpeg("CountyRisk_overMin_hist.jpeg", width = 760, height = 520, units = 'px', res = 100)
    hist(CountyRisk.long$Value[CountyRisk.long$Value > 0], breaks = 25, xlab = "County Risk",
         main = paste0("County Risk > ",CountyRisk_min))
    dev.off()
    
    
    # violin plot
    jpeg("CountyRisk_overMin_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = CountyRisk.long.overMin, aes(x = CountyRisk.long.overMin$type, y = CountyRisk.long.overMin$Value, color = CountyRisk.long.overMin$type)) +
        scale_x_discrete(limits = levels(CountyRisk.long.overMin$type)) + theme_bw() +
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        labs(x = NULL, y = "County Risk", title = paste0("County Risk > ",CountyRisk_min)) + scale_color_manual(values = cbPalette)
    }) 
    dev.off()
    
    jpeg("CountyRisk_Low_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = CountyRisk.long.low, aes(x = CountyRisk.long.low$type, y = CountyRisk.long.low$Value, color = CountyRisk.long.low$type)) +
        scale_x_discrete(limits = levels(CountyRisk.long.low$type)) + theme_bw() +
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        scale_color_manual(values = cbPalette) + labs(x = NULL, y = "County Risk",
                                                      title = paste0("County Risk \u2264 ",as.character(CountyRisk_cutoff)))
    })
    dev.off()
    
    jpeg("CountyRisk_High_violin.jpeg", width = 760, height = 520, units = 'px', res = 100)
    print({
      ggplot() + geom_violin(data = CountyRisk.long.high, aes(x = CountyRisk.long.high$type, y = CountyRisk.long.high$Value, color = CountyRisk.long.high$type)) +
        scale_x_discrete(limits = levels(CountyRisk.long.high$type)) + theme_bw() +
        theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        scale_color_manual(values = cbPalette) + labs(x = NULL, y = "County Risk",
                                                      title = paste0("County Risk > ",as.character(CountyRisk_cutoff)))
    })
    dev.off()
    
    
    CountyRisk_Median_RepType=tapply(CountyRisk.long$Value,CountyRisk.long$type,FUN=median)
    CountyRisk_Upper_RepType<-tapply(CountyRisk.long$Value,CountyRisk.long$type,FUN=quantile,probs=0.975, na.rm = TRUE)
    
    
    
    ## Create individual data frames for each run type. Add medians and upper 2.5% quartiles, then split into
    # two-column dataframes with fips for mapping
    
    CountyRisk_medians=NULL
    CountyRisk_uppers=NULL
    
    for (i in 1:length(run.types)){
      
      df <- CountyRisk[,c(1,2,((i-1)*2*runs_per_ctrl_type+3):(i*2*runs_per_ctrl_type+2))]
      df1 <- df[df[,4] != 0 & !is.na(df[,4]),]
      newname <- df1[1,4]
      df.new <- df[,-grep("type_", colnames(df))]
      
      # Calculate median and upper 2.5%
      df.new$median=apply(df.new[grep("run_", colnames(df.new))], 1, median, na.rm = TRUE)
      df.new$upper=apply(df.new[grep("run_", colnames(df.new))], 1, quantile, probs=0.975, na.rm = TRUE)
      
      # Create the vectors for map scales
      CountyRisk_medians=c(CountyRisk_medians,df.new$median)
      CountyRisk_uppers=c(CountyRisk_uppers,df.new$upper)
      
      # save datasets for mapping, one for each level with medians and uppers with name of control type
      assign(paste0("CountyRisk_med_",newname), df.new[,grepl( "fips|median" , names( df.new ) )])
      assign(paste0( "CountyRisk_upper_",newname), df.new[,grepl( "fips|upper" , names( df.new ) )])
      
      rm(df, df1,df.new,newname)
    }
    
    ## Median maps
    
    # Create a vector for scale for median maps
    CountyRisk_median_values <- round(c(0, 0.001, mean(CountyRisk_medians, na.rm = TRUE),
                                        mean(CountyRisk_medians, na.rm = TRUE) + sd(CountyRisk_medians, na.rm = TRUE),
                                        mean(CountyRisk_medians, na.rm = TRUE) + 3*sd(CountyRisk_medians, na.rm = TRUE),
                                        max(CountyRisk_medians, na.rm = TRUE)), 3)
    CountyRisk_median_values<-sort(unique(CountyRisk_median_values))
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("CountyRisk_med_",run.types[i]))
      jpeg(paste0(path_output, paste0("CountyRisk_Median_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30",
                  missing.include = TRUE, color.break.type = "values", legend.digits = 3,
                  color.break.values = CountyRisk_median_values, color.sequence = palette,
                  legend.spacing = 6.5, legend.shrink = 0.3, legend.width = 1)
      dev.off()
    }
    
    ## Upper 2.5% maps
    
    # Create a vector for scale for uppers maps
    CountyRisk_upper_values <- round(c(0, 0.001, mean(CountyRisk_uppers, na.rm = TRUE),
                                       mean(CountyRisk_uppers, na.rm = TRUE) + sd(CountyRisk_uppers, na.rm = TRUE),
                                       mean(CountyRisk_uppers, na.rm = TRUE) + 3*sd(CountyRisk_uppers, na.rm = TRUE),
                                       max(CountyRisk_uppers, na.rm = TRUE)), 3)
    CountyRisk_upper_values<-sort(CountyRisk_upper_values)
    
    for (i in 1:length(run.types)){
      name_df=get(paste0("CountyRisk_upper_",run.types[i]))
      jpeg(paste0(path_output, paste0("CountyRisk_Up_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30",
                  missing.include = TRUE, color.break.type = "values", legend.digits = 3,
                  color.break.values = CountyRisk_upper_values, color.sequence = palette,
                  legend.spacing = 6.5, legend.shrink = 0.3, legend.width = 1)
      dev.off()
    }
  }
  
  #####################################################################################################
  #####################################################################################################
  
  #####  Type of spread aka proportion local transmission (Detail Files)  #######
  
  #####################################################################################################
  #####################################################################################################
  
  if (localSpread == TRUE){
    
    # Extract the infection route from the detail file
    setwd(pathfiles)
    
    detail.res <- as.data.frame(fread(detail.fnames[1], header = TRUE, select = c("SourceCounty", "InfRoute", "ControlPrevented")))
    detail.res <- detail.res[detail.res$ControlPrevented == "none",]
    detail.res <- as.data.frame(table(detail.res$SourceCounty, detail.res$InfRoute))
    TypeSpread <- dcast(as.data.table(detail.res), Var1~Var2, value.var = "Freq")
    TypeSpread <- as.data.frame(TypeSpread)
    
    colnames(TypeSpread) <- c("SourceCounty", "Local", "Ship")
    
    TypeSpread$Type <- unlist(strsplit(detail.fnames[1], "_FLAPS"))[1]

    for(i in 2:length(detail.fnames)){
      detail.res <- as.data.frame(fread(detail.fnames[i], header = TRUE, select = c("SourceCounty", "InfRoute", "ControlPrevented")))
      detail.res <- detail.res[detail.res$ControlPrevented == "none",]
      detail.res <- as.data.frame(table(detail.res$SourceCounty, detail.res$InfRoute))
      detail.res <- dcast(as.data.table(detail.res), Var1~Var2, value.var = "Freq")
      detail.res <- as.data.frame(detail.res)
      
      colnames(detail.res) <- c("SourceCounty", "Local", "Ship")
      
      detail.res$Type <- unlist(strsplit(detail.fnames[i], "_FLAPS"))[1]
      TypeSpread <- merge(TypeSpread, detail.res, by = "SourceCounty", all=TRUE)
      
    }
    
    names(TypeSpread)[seq(2,ncol(TypeSpread),3)] <- paste0("local_", sprintf("%03d", seq(1:((ncol(TypeSpread))/3))))
    names(TypeSpread)[seq(3,ncol(TypeSpread),3)] <- paste0("ship_", sprintf("%03d", seq(1:((ncol(TypeSpread))/3))))
    names(TypeSpread)[seq(4,ncol(TypeSpread),3)] <- paste0("type_", sprintf("%03d", seq(1:((ncol(TypeSpread))/3))))
    
    setwd(path_output)
    
    if (export.datafiles == 1 | export.datafiles == 3) {write.csv(TypeSpread, "TypeSpread.csv")}

    
    ## Create individual data frames for each run type. Add proportions local and shipment spread, then split into 
    # a two-column dataframe with fips and proportion local spread for mapping.  
    
    for (i in 1:length(run.types)){
      
      df <- TypeSpread[,c(1,((i-1)*3*runs_per_ctrl_type+2):(i*3*runs_per_ctrl_type+1))]
      newname <- first(na.omit(df[,4]))
      df.new <- df[,-grep("type_", colnames(df))]
      
      # Calculate prop local and shipment 
      df.new$local=apply(df.new[grep("local_", colnames(df.new))], 1, sum, na.rm = TRUE)
      df.new$ship=apply(df.new[grep("ship_", colnames(df.new))], 1, sum, na.rm = TRUE)
      
      # Calculate proportion local 
      df.new$propLocal <- round((df.new$local/(df.new$local + df.new$ship))*100,0)
      
      # Generate two-column df for mapping
      assign(paste0("PropLocal_",newname), df.new[,grepl( "SourceCounty|propLocal" , names( df.new ) )])
      
      # Cleanup
      rm(df, df.new,newname)
    }
    
    ## Maps of the proportion of spread events that are local spread. 
    # Create a vector for scale
    local_scale <- round(c(0, 25, 50, 60, 70, 80, 85, 90, 95, 97.5, 99, 100),2)
    
    # loop over run types
    for (i in 1:length(run.types)){
      name_df=get(paste0("PropLocal_",run.types[i]))
      jpeg(paste0(path_output, paste0("PropLocal_",run.types[i],".jpeg"), sep=""), width = 760, height = 520, units = 'px', res = 100)
      setwd(path0)
      
      map_by_fips(name_df, county.border.col = NA, state.border.col = "gray30",
                  missing.include = TRUE, color.break.type = "values", legend.digits = 1,
                  color.break.values = local_scale, 
                  color.sequence = if(ls_match == TRUE) {palette} else {color_bluepurple}, 
                  legend.spacing = 5.5, legend.shrink = 0.6, legend.width = 1)
      dev.off()
    }
    
  }

  
  #####################################################################################################
  #####################################################################################################
  
  ####### Summary table  #######
  
  #####################################################################################################
  #####################################################################################################
  
  if (summaryTable == TRUE){
    
    num.metrics = sum(exists("Anim_Median_RepType"),exists("Dur_Median_RepType"),exists("EpidExt_Median_RepType"),exists("PremInf_Median_RepType"))
    
    metrics = NULL
    if(exists("Dur_Median_RepType")) {metrics <- c(metrics,"Dur_Median_RepType","Dur_Upper_RepType")}
    if(exists("PremInf_Median_RepType")) {metrics <- c(metrics,"PremInf_Median_RepType","PremInf_Upper_RepType")}
    if(exists("Anim_Median_RepType")) {metrics <- c(metrics,"Anim_Median_RepType","Anim_Upper_RepType")}
    if(exists("EpidExt_Median_RepType")) {metrics <- c(metrics,"EpidExt_Median_RepType","EpidExt_Upper_RepType")}
    
    metric.alias = NULL
    if(exists("Dur_Median_RepType")) {metric.alias <- c(metric.alias,"Duration")}
    if(exists("PremInf_Median_RepType")) {metric.alias <- c(metric.alias,"Premises Infected")}
    if(exists("Anim_Median_RepType")) {metric.alias <- c(metric.alias,"Infected Animals")}
    if(exists("EpidExt_Median_RepType")) {metric.alias <- c(metric.alias,"Epidemic Extent (counties)")}
    
    sumtab=NULL
    
    for(i in 1:length(metrics)){
      m <-get(metrics[i])
      sumtab=cbind(sumtab,m) 
    }
    
    colnames(sumtab) <- rep(c("Median", "Upper 2.5%"),num.metrics)
    rownames(sumtab) <- run.types
    
    aboveheader = c( 1 , rep(2, times = length(metric.alias)))
    names(aboveheader ) = c( "" , metric.alias )
    
    setwd(path_output)
    save_kable(x= kable(sumtab, "html") %>%
                 kable_styling(bootstrap_options = c("striped", "hover", "condensed", "bordered"), full_width = F, position = "float_right") %>%
                 add_header_above(aboveheader),
               file=paste0("Summary_Table_",format(Sys.time(),'%Y%m%d_%H%M'),".jpeg"))
    
    
  }
  
}




processUSDOS(summaryTable = T,
             duration = TRUE,
             premInf = F,
             premReport = F,
             epidemicExtent = F,
             movementBan = F,
             premisesCulled = F,
             premisesVax = F,
             animalsInfected = F,
             countyRisk = F,
             localSpread = F,
             color_palette = "color_red",
             ls_match = TRUE)



