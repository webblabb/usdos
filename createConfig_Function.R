
#### Pre-Processing Function for USDOSv2


#####################################################################################################################################
## Working directory must contain "templates" folder containing the SH_FILE_TEMPLATE.txt as well as the template to be used for the 
##   config files. It also must either be or contain a destination folder for the config files that will be generated. We call this folder 
##   "Generated_configs".

#####################################################################################################################################
## In addition to config, batch, and job files, this code generates a dataframe, saved as a .csv, with the following information:
##   The first column contains the filenames of the configs that will be generated
##   The second column contains the jobnames that will be displayed as each file is running
##   All other columns contain the information to be inserted into the config files, with the column headers matching the
##     appropriate line number in the config file.

#####################################################################################################################################
#' USDOS: Generate config files
#' 
#' A streamlined method for producing large numbers of USDOS config files and changing parameters.
#' 
#' @param sim_type The initial parameter to assign the control to be initiated. This value changes the default values for specific additional parameters. The default value is 'base', indicating no controls will be implemented and the default values are as listed. Other options include 'IPcull',IPDCcull', 'cullVax', 'vax', and 'sensitivity'. 'IPcull' will by default change the 'ctrl_type' to 1, indicating cull of infected premises. 'IPDCcull' will by default change the 'ctrl_type' to 2, indicating cull of infected premises and dangerous contacts, and will set vax_range to -1 to indicate an action for DCs.Selecting 'cullVax' will change the 'ctrl_type' to 3, 'vax_eff' to 1, 'vax_delay' to 11, and 'vax_range' to NULL. Selecting 'vax' will change the 'ctrl_type' to 4, 'vax_eff' to 1, 'vax_delay' to 11, and 'vax_range' to NULL. 'sensitivity' will change 'source' to xx and 'filetype' to xx.
#' @param template.config.location The location of the template for the config file. Default is 'templates/USDOS_CONFIG_TEMPLATE_NoControls_workaround.txt', a folder within the current working directory. This will be automatically changed with each 'sim_type'.
#' @param sh_template The destination of the SH_FILE_TEMPLATE.txt file. Default is 'templates/SH_FILE_TEMPLATE.txt'.
#' @param jobfile.name The name of the jobfile to output. Default is 'USDOS_Sim.job'. 
#' @param reps The number of times each FLAPS will be repeated. The default value is 10, for a total of 100 runs.
#' @param cutoff The max number of infection events before stopping each run. '*' indicates no limit to the infection
#' @param cutoff_days The number of timesteps (days) to run for each seeded outbreak. Default is 365.
#' @param source The source FIPS to be seeded, filename, or 'allFips'. Filename is of the tab-delimited file containing the FIPS codes or premisesIDs from which to seed infection, one line per simulation. 'allFips' will seed from all counties containing premises. Seeding from multiple premises at once requires comma-separated premisesIDs in a file, one line per simulation. Default is 'allFips'.
#' @param filetype The seed file type - type of information provided per line of source. 'fips' will choose 1 premises at random within the FIPS codes listed in the source parameter or in all FIPS, 'singlePremises' will use the premisesID provided by the source parameter, 'multiplePremises' will seed all comma separated premisesIDs in the source parameter. Default is 'fips'.
#' @param parameter.sample If "spread" or "control", a set of parameters to be used in analysis is provided in a supplemental file called "LHS_Parameters.csv". "spread" is used to provide variation in Tc, K2, K3, latency, and infectious parameters. "control" is used to provide variation in cull_rate, vax_rate, vax_delay, index_rep_time, rep_time, DC_rep_time, DC_sus, and DC_exp parameters. The supplemental file must contain a number of rows equal to the total number of simulations being generated. Default is NA.
#' @param fips.info Name of file containing fips name, state name, area (m2), x, y. Tab separated. Default is "inputfiles/FIPS_20151805.txt"
#' @param flaps The subname of the FLAPS version being used. For example, '12_min' will use the files 'flaps12_min_0001.txt' to 'flaps12_min_0010.txt'.
#' @param landfill.files Set to TRUE if there is a unique landfill file for each FLAPS file. TRUE will require landfill files named "landfill_{filename of flaps}.txt" (e.g. 'landfills_flaps12_min_0001.txt'. Otherwise, set to FALSE and the file name listed in the config template file will be used. If the movement ban has a file associated, set to false and manually enter into the config file template.
#' @param vaccine.file Set to the filename if there is a vaccination constraint file. The default for when vaccination is on will be a vaccination file named "vaccineBankUpdated.txt" in the inputfiles folder. Otherwise, set to "NA" for no vaccination file.
#' @param shipment.gen Method to generate county-to-county shipments with USAMM (1 = USAMM version 1, -1 = shipments off)
#' @param Tc The infectiousness constant. Default is 10.252.
#' @param k1 Kernel parameter 1. Default is 1.46e-08.
#' @param k2 Kernel parameter 2. Default is 1686.155.
#' @param k3 Kernel parameter 3. Default is 2.267.
#' @param susc.exp.sp1 Susceptibility exponent (power, mu) for the first species listed on line 12 of the config file. Default species is beef, with power = 1 
#' @param susc.exp.sp2 Susceptibility exponent (power, mu) for the second species listed on line 12 of the config file. Default species is dairy, with power = 1 
#' @param infect.exp.sp1 Infectiousness exponent (power, omega) for the first species listed on line 12 of the config file. Default species is beef, with power = 1 
#' @param infect.exp.sp2 Infectiousness exponent (power, omega) for the first species listed on line 12 of the config file. Default species is dairy, with power = 1 
#' @param latency The mean days from premises exposure to infectiousness. Default is 5 days.
#' @param infectious The mean days from premises infectiousness to immunity. Default is 7 days.
#' @param ban_eff The effectiveness of the movement ban put into place on shipments. The value is between 0 and 1, with a default value of 0 for no movement ban. 1 indicates a 100% effective movement ban.
#' @param ctrl_type The type of control to be implemented in the runs. Default is 0 for no control. Other values are 1 for infected premises cull, 2 for infected premises cull and dangerous contact cull, 3 for infected premises cull and either dangerous contact or ring vaccination, and 4 for infected premises vaccination.
#' @param vax_eff The effectiveness of the vaccination. The value is between 0 and 1, and is only in effect if a vaccination control type is selected. The default value is 1 for 100% effectiveness.
#' @param vax_delay The length of time, in days, between vaccination and vaccine effectiveness. This is only in effect if a vaccination control type is selected.
#' @param vax_range The range of the vaccination ring. The default value is '-1', indicating that the vaccination is applied to dangerous contacts, and not in a ring. A positive number indicates a ring vaccination, in meters. This is only in effect if a vaccination control type is selected.
#' @param cull_rate The rate at which species are culled on a premises. The default is 240 cattle per day. 
#' @param vax_rate The rate at which species are vaccinated on a premises. The default is 6804 cattle per day.
#' @param DC_sus The scale for tracing DCs of a susceptible premises. The default value is 4 to make the mean DCs per Reported premises between 1.5-2 when combined with the exposed value.
#' @param DC_exp The scale for tracing DCs of a exposed premises. The default value is 5 to make the mean DCs per Reported premises between 1.5-2 when combined with the susceptible value.
#' @param index_rep_time The mean time (days) from index premises exposure to reporting. Default is 15 days.
#' @param rep_time The mean time (days) from non-dangerous contact premises exposure to reporting. Default is 8 days.
#' @param DC_rep_time The mean time (days) from dangerous contact premises exposure to reporting. Default is 2 days.
#' @param destination.folder The location to put all generated files. Folder must exist. Default is the current working directory.
#' @param replacement.df The dataframe containing all information to use to modify the default config file. This is created with the time, date, and replacement_df as its title.
#' @param sens.number The row number in the parameter set document that you want to use for sensitivty.
#'
#' @return This code will generate a config file and batch file for each of the lines in the replacement.df dataframe.
#' @author Send bug reports, suggestions, corrections, or comments to Deedra Murrieta or Katie Owers (kowers@colostate.edu).

#' @export

#

createConfigs <- function(sim_type = "base", ...){
  
  ### Set default values for config files:
  
  # Inputs
  template.config.location <- "templates/USDOS_CONFIG_TEMPLATE_NoControls_workaround.txt"
  sh_template <- "templates/SH_FILE_TEMPLATE.txt"
  jobfile.name <- "USDOS_Sim.job"
  
  # Simulation run time
  reps <- 10
  cutoff <- "*"
  cutoff_days<- 365
  flaps <- "12_min"
  landfill.files <- TRUE
  vaccine.file <- "NA"
  shipment.gen <- 1
  fips.info <- "inputfiles/FIPS_20151805.txt"
  
  # Sensitivity 
  source <- "allFips" 
  filetype <- "fips"
  parameter.sample <- "NA"
  sens.number <- NA
  
  # Infection Parameters
  Tc <- 10.252
  k1 <- 1.46e-08
  k2 <- 1686.155
  k3 <- 2.267
  latency <- 5
  infectious <- 7
  susc.exp.sp1 <- 1
  susc.exp.sp2 <- 1
  infect.exp.sp1 <- 1
  infect.exp.sp2 <- 1
  
  # Control Parameters
  ban_eff <- 0
  ctrl_type <- 0
  cull_rate <- 240
  vax_rate <- NULL
  vax_eff <- 0
  vax_delay <- NULL 
  vax_range <- NULL
  
  # DC Reporting
  index_rep_time <- 15
  rep_time <- 8
  DC_rep_time <- 2
  DC_sus <- 4
  DC_exp <- 5
  
  # Outputs
  destination.folder <- getwd()
  replacement.df <- data.frame()
  
  
  ### Set the default simulation types with if statements.
  if(sim_type == "IPcull"){
    ## Change default settings for IP cull runs
    ctrl_type = 1 
    template.config.location <- "templates/USDOS_CONFIG_TEMPLATE_MovementBan_IP.txt"
    
  }
  
  if(sim_type == "IPDCcull"){
    ## Change default settings for IP + DC cull runs
    ctrl_type = 2 
    vax_range = -1
    template.config.location <- "templates/USDOS_CONFIG_TEMPLATE_MovementBan_IP_DC.txt"
    
  }
  
  if(sim_type == "cullVax"){
    ## Change default settings for cull & vaccination runs
    ctrl_type <- 3
    vax_rate <- 6804
    vax_eff <- 0.9
    vax_delay <- 11 
    vax_range <- -1 # -1 corresponding to DC vax
    vaccine.file <- "vaccineBankUpdated.txt"
    template.config.location <- "templates/USDOS_CONFIG_TEMPLATE_MovementBan_IP_VAX.txt"
    
  }
  
  # if(sim_type == "vax"){
  #   ## Change default settings for vaccination runs
  #   ctrl_type = 4
  #   vax_eff = 1
  
  #   vax_delay = 11 
  #   vax_range = NULL
  #   template.config.location = "templates/USDOS_CONFIG_TEMPLATE_MovementBan_VAX.txt"
  # } 
  
  if(sim_type == "sensitivity"){
    ## Change default settings for sensitivity runs
    source <- "inputfiles/SensCounty_shortlist_0219.txt"
    filetype <- "fips"
  }
  
  parameters = c("reps", "cutoff", "cutoff_days", "source", "filetype", "parameter.sample", "shipment.gen",
                 "fips.info", "Tc", "k1", "k2", "k3", "latency", "infectious", "ban_eff", "cull_rate", "vax_rate",
                 "vax_eff", "vax_delay", "vax_range", "index_rep_time", "rep_time", 
                 "DC_rep_time", "DC_sus", "DC_exp","sens.number", "jobfile.name", "landfill.files", "vaccination.file")
  
  ## Override values based on additional arguments passed to the function. 
  ## Add names to this list if they are additional arguments for the function.
  additional.inputs = list(...)
  parameter.options = c("template.config.location", "sh_template", "jobfile.name", "flaps", "landfill.files", "vaccination.file",
                        "fips.info", "reps", "cutoff", "cutoff_days", "source", "filetype", "parameter.sample", 
                        "shipment.gen", "Tc", "k1", "k2", "k3","susc.exp.sp1","susc.exp.sp2","infect.exp.sp1" ,
                        "infect.exp.sp2","latency", "infectious", "ban_eff", "cull_rate", "vax_rate", "ctrl_type", 
                        "vax_eff", "vax_delay", "vax_range", "index_rep_time", "rep_time", 
                        "DC_rep_time", "DC_sus", "DC_exp","destination.folder", "replacement.df", "sens.number")
  
  ## Return error if additional inputs aren't potential parameters
  unknown.inputs = additional.inputs[!names(additional.inputs) %in% parameter.options]
  if(length(unknown.inputs) > 0){
    error.message = paste0("----------------------------------------------------\n",
                           "These arguments aren't recognized: \n",
                           paste(names(unknown.inputs), collapse = "\n"), 
                           "\n",
                           "----------------------------------------------------\n")
    return(cat(error.message))
  }
  ## Use only the additional arguments that match the possible parameter names
  additional.inputs = additional.inputs[names(additional.inputs) %in% parameter.options]
  
  for(parameter in names(additional.inputs)){
    assign(parameter, additional.inputs[parameter][[1]]) 
  }
  
  ## Create different parameter values based on the latin hypercube sampling
  # The csv with the parameter list needs to be in the inputfiles folder.
  if(parameter.sample == "spread"){
    lh <- read.csv("inputfiles/LHS_Parameters.csv", header = TRUE)
    Tc <- lh$Tc[sens.number]
    k1 <- lh$k1[sens.number]
    k2 <- lh$k2[sens.number]
    k3 <- lh$k3[sens.number]
    latency <- lh$lat[sens.number]
    infectious <- lh$inf[sens.number]
  }
  if(parameter.sample == "control"){
    lh <- read.csv("inputfiles/LHS_Control_Parameters_040119.csv", header = TRUE)
    cull_rate <- lh$cull_rate[sens.number]
    vax_rate <- lh$vax_rate[sens.number]
    vax_delay <- lh$vax_delay[sens.number]
    index_rep_time <- lh$index_rep_time[sens.number]
    rep_time <- lh$rep_time[sens.number]
    DC_rep_time <- lh$DC_rep_time[sens.number]
    DC_sus <- lh$DC_sus[sens.number]
    DC_exp <- lh$DC_exp[sens.number]
    ban_eff <- lh$shipban[sens.number]
    if(sim_type == "IPDCcull" | sim_type == "IPcull"){vax_rate = NULL}
    if(sim_type == "IPDCcull" | sim_type == "IPcull"){vax_delay = NULL}
  } 
  
  ## Create list of FLAPS files from the inputted flaps argument. 
  ## FLAPS name is an argument input above, default is "12", we are likely to use "12_min"
  flaps.names <- paste0("flaps",flaps,"_", sprintf("%04d", 1:10), ".txt", sep = "")
  all.flaps.names <- rep(flaps.names, each = reps)
  
  ## Create list of additional control files from the inputted control argument  
  if(landfill.files == TRUE){
    landfill.file.list <- paste0("landfills_flaps",flaps,"_", sprintf("%04d", 1:10), ".txt", sep = "")
    landfill.file.names <- rep(landfill.file.list, each = reps)}
  
  ## Create file naming conventions based on: 
  
  # control type selections
  type <- if(ctrl_type == 0){"base_"}else if(ctrl_type == 1){"IP_"}else if(ctrl_type == 2){"IP_DC_"}else if(ctrl_type == 3){"IP_VAX_"}else if(ctrl_type == 4){"VAX_"}else{"_"}
  
  # vaccine ring radius (also includes check for allowable values of vax_range)
  vax_rng <-if(is.null(vax_range) == T){-1}else if(vax_range == -1){-1} else if(vax_range>=0){vax_range} else stop ("Vaccination Range must be positive or '-1' if not using ring vaccination")
  
  vax <- if(vax_rng>=0){paste0(vax_rng/1000,"km_")}else{""} 
  
  # movement ban effectiveness
  MvmtBan <- if(ban_eff!=0){paste0("MvmtBan_",100*ban_eff,"_")}else{""}
  
  # Shipments off
  ShipOff <-if(shipment.gen == -1) {"ShipmentsOff_"} else{""}
  
  # Errors for parameters that are out of bounds and warnings for unique situations
  if(ban_eff < 0 | ban_eff > 1) stop ('Movement Ban Effectiveness must be between 0 and 1')
  if(vax_eff < 0 | vax_eff > 1) stop ('Vaccine Effectiveness must be between 0 and 1')
  
  
  ## Create the replacement_df data frame that contains all of the information created by this function. 
  ## There needs to be the ability to add columns independently when things are manually changed in the function.
  ## We will need to add the variables that are added automatically with the control type and the manual ones.
  config.fname <- paste0("config_", 
                         type,vax,MvmtBan,ShipOff,
                         gsub(".txt", "", all.flaps.names), "_", sprintf("%02d", rep(1:10, 10)), 
                         format(Sys.time(), '_%Y%m%d'),"_", ifelse(is.na(sens.number), "", sens.number), ".txt")
  
  jobname <- paste0("MI_",type,vax,MvmtBan, ShipOff,"f", sprintf("%02d", rep(1:10, each=10)), "_", 
                    sprintf("%02d", rep(1:10), "_", sens.number))
  
  replacement.df <- data.frame(config.fname, jobname, 
                               FLAPS.locations = paste0("FLAPS/", all.flaps.names), 
                               resources = if(landfill.files == TRUE & vaccine.file == "NA" & vax_eff == 0){resources = paste0("NA;inputfiles/", landfill.file.names)}else 
                                 if(landfill.files == TRUE & vaccine.file != "NA"){resources = paste0("NA;inputfiles/", landfill.file.names, ";inputfiles/", vaccine.file)}else
                                    if(landfill.files == FALSE & vaccine.file != "NA"){resources = paste0("NA;NA;inputfiles/", vaccine.file)}else
                                   {paste0(strsplit(readLines(template.config.location)[74], "#")[[1]][1])},
                               resource.filetype = if(landfill.files == TRUE & vaccine.file == "NA" & vax_eff == 0){resource.filetype = paste0("NA;resourceLocs")}else 
                                 if(landfill.files == TRUE & vaccine.file != "NA"){resources = paste0("NA;resourceLocs;resourceBoosts")}else
                                   if(landfill.files == FALSE & vaccine.file != "NA"){resources = paste0("NA;NA;resourceBoosts")}else
                                   {paste0(strsplit(readLines(template.config.location)[75], "#")[[1]][1])},
                               cutoff = paste0(cutoff),
                               cutoff_days = paste0(cutoff_days),
                               source = paste0(source), 
                               filetype = paste0(filetype),
                               fips.info = paste0(fips.info),
                               shipment.gen = paste0(shipment.gen),
                               Tc = paste0(Tc, ",", Tc), 
                               kernels = paste0(k1,",", k2, ",", k3), 
                               susceptibility.exponents=paste0(susc.exp.sp1,",",susc.exp.sp2),
                               infectiousness.exponents=paste0(infect.exp.sp1,",",infect.exp.sp2),
                               latency = paste0(latency, ",0"), 
                               infectious = paste0(infectious, ",0"),
                               shipban.effectiveness = if(vax_eff == 0){paste0(ban_eff, ",", ban_eff, ";1,1")}else{
                                 paste0(ban_eff, ",", ban_eff, ";1,1;", vax_eff, ",", vax_eff)}, 
                               constraint = if(is.null(vax_rate)){paste0("0;", cull_rate, ",0")}else{
                                 paste0("0;", cull_rate, ",0;", vax_rate, ",0")},
                               vax_delay = if(is.null(vax_delay)){paste0("0,0")}else{paste0("0,0,", vax_delay)}, 
                               vax_range = if(is.null(vax_range)){paste0("0,0")}else{paste0("0,0,", vax_rng)}, 
                               index_rep_time = paste0(index_rep_time, ",0"),
                               rep_time = paste0(rep_time, ",0"),
                               DC_rep_time = paste0(DC_rep_time, ",0"),
                               DC_scaling = paste0("sus,", DC_sus, ";exp,", DC_exp),
                               stringsAsFactors = FALSE)
  
  
  replacement.df$batch.name = gsub("config_", "", replacement.df$config.fname)
  replacement.df$batch.name = gsub(".txt", "", replacement.df$batch.name)
  
  columns <- length(replacement.df)
  for(i in 3:columns){
    if (colnames(replacement.df)[i] == "batch.name"){colnames(replacement.df)[i] = 1} else 
      if(colnames(replacement.df)[i] == "FLAPS.locations"){colnames(replacement.df)[i] = 11} else 
        if (colnames(replacement.df)[i] == "cutoff_days"){colnames(replacement.df)[i] = 13} else 
          if (colnames(replacement.df)[i] == "cutoff"){colnames(replacement.df)[i] = 14} else 
            if (colnames(replacement.df)[i] == "fips.info"){colnames(replacement.df)[i] = 18} else
              if (colnames(replacement.df)[i] == "source"){colnames(replacement.df)[i] = 21} else 
                if (colnames(replacement.df)[i] == "filetype"){colnames(replacement.df)[i] = 22} else 
                  if (colnames(replacement.df)[i] == "susceptibility.exponents"){colnames(replacement.df)[i] = 24} else
                    if (colnames(replacement.df)[i] == "infectiousness.exponents"){colnames(replacement.df)[i] = 25} else
                      if (colnames(replacement.df)[i] == "Tc"){colnames(replacement.df)[i] = 27} else
                        if (colnames(replacement.df)[i] == "kernels"){colnames(replacement.df)[i] = 29} else
                          if (colnames(replacement.df)[i] == "latency"){colnames(replacement.df)[i] = 31} else
                            if (colnames(replacement.df)[i] == "infectious"){colnames(replacement.df)[i] = 32} else
                              if (colnames(replacement.df)[i] == "shipment.gen"){colnames(replacement.df)[i] = 41} else 
                                if (colnames(replacement.df)[i] == "constraint"){colnames(replacement.df)[i] = 53} else
                                  if(colnames(replacement.df)[i] == "resources"){colnames(replacement.df)[i] = 55} else 
                                    if(colnames(replacement.df)[i] == "resource.filetype"){colnames(replacement.df)[i] = 56} else
                                      if (colnames(replacement.df)[i] == "vax_delay"){colnames(replacement.df)[i] = 57} else
                                        if (colnames(replacement.df)[i] == "shipban.effectiveness"){colnames(replacement.df)[i] = 61} else
                                          if (colnames(replacement.df)[i] == "vax_range"){colnames(replacement.df)[i] = 66} else
                                            if (colnames(replacement.df)[i] == "index_rep_time"){colnames(replacement.df)[i] = 71} else
                                              if (colnames(replacement.df)[i] == "rep_time"){colnames(replacement.df)[i] = 72} else
                                                if (colnames(replacement.df)[i] ==  "DC_rep_time"){colnames(replacement.df)[i] = 73} else
                                                  if (colnames(replacement.df)[i] == "DC_scaling"){colnames(replacement.df)[i] = 74} else
                                                    warning(paste0(colnames(replacement.df)[i]," has not been added to this function. Change the line 
                                                       value in the template file then rerun this function without that variable."))
  }
  
  ## Write out the data frame with the associated parameters
  print(head(replacement.df))
  str1 <- 'replacementDF.csv' 
  write.csv(replacement.df, paste0(format(Sys.time(),'%Y%m%d_%H%M_'), sub('\\..*', '', str1), ".csv"), row.names = FALSE)
  #write.csv(replacement.df, paste0(destination.folder,"/",paste0(format(Sys.time(),'%Y%m%d_%H%M_'), sub('\\..*', '', str1), ".csv")), row.names = FALSE) # For moving to subdirectory in Windows
  
  ## Generate the Config, Batch, and Job files.
  destination.folder = paste0(destination.folder, "/")
  config = readLines(template.config.location)
  sh.template = readLines(sh_template)
  for(i in 1:dim(replacement.df)[1]){  
    config = readLines(template.config.location)
    for(j in 3:dim(replacement.df)[2]){
      replacement.value = replacement.df[i,j]  
      # Find the appropriate line:
      replacement.expression = paste0("\\#\\(", colnames(replacement.df)[j], ")")
      replacement.line.no = grep(replacement.expression, config)
      
      replacement.line = config[replacement.line.no]
      comment.string = strsplit(replacement.line, replacement.expression)[[1]][2]
      config[replacement.line.no] = 
        paste0(replacement.value," ", gsub("\\\\", "", replacement.expression), comment.string)
      
    }
    fname = replacement.df$config.fname[i]
    cat(config, file=paste0(destination.folder, fname), sep="\n")
    
    jobname = replacement.df$jobname[i]
    
    #Generate a batch (.sh) file for each config:
    sh.file = gsub("JOBNAME", jobname, sh.template)
    sh.file = gsub("CONFIGNAME", fname, sh.file)
    fname.sh = gsub(".txt", ".sh", fname)
    fname.sh = gsub("config", "BATCH",fname.sh)
    cat(sh.file, file = paste0(destination.folder, fname.sh), sep = "\n")
    
    ### Add a line to the jobfile so the Batch file will run:
    cat(paste0("sbatch ", fname.sh, "\n"), file=paste0(destination.folder, jobfile.name), append=TRUE)
  }
}


## Example commands

#### No Control (base run) ####
# createConfigs(sim_type = "base", jobfile.name="base.job",
#             flaps = "12_min",landfill.files = TRUE)

##### IP cull, 90% effective movement ban ####
# createConfigs(sim_type = "IPcull", ban_eff = 0.9, jobfile.name="IP_90.job", 
#              flaps = "12_min", landfill.files = TRUE)

#### IP + DC cull, 90% effective movement ban ####
# createConfigs(sim_type = "IPDCcull", ban_eff = 0.9, jobfile.name="IP_DC_90.job",
#               flaps = "12_min", landfill.files = TRUE)

#### IP cull + DC vax, 90% effective movement ban ####
# createConfigs(sim_type = "cullVax",ban_eff = 0.9, jobfile.name="IP_VAX_90.job",
#              flaps = "12_min", landfill.files = TRUE)

#### IP cull + 3km ring vax, 90% effective  movement ban ####
# createConfigs(sim_type = "cullVax", vax_range = 3000, ban_eff = 0.9, jobfile.name="IP_VAX_90_3.job",
#              flaps = "12_min", landfill.files = TRUE)

#### IP cull + 10km ring vax, 90% effective movement ban ####
# createConfigs(sim_type = "cullVax", vax_range = 10000, ban_eff = 0.9, jobfile.name="IP_VAX_90_10.job",
#              flaps = "12_min", landfill.files = TRUE)

# # Sensitivity
# for(i in 1:100){
# 
# createConfigs(source = "inputfiles/SensCounty_shortlist_0219.txt", filetype = "fips", sim_type = "IPDCcull", 
#               parameter.sample = "control", flaps = "12_min",  sens.number = i, 
#               jobfile.name = paste0("USDOS_Sens_", i, ".job"))
# }




