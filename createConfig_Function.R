#### USDOS Pre-Processing Function 

# Code by Deedra Murrieta & Katie Owers,

#####################################################################################################################################
## Working directory must contain "templates" folder containing the SH_FILE_TEMPLATE.txt as well as the template to be used for the 
##   config files. It also must either be or contain a destination folder for the files that will be generated. 

#####################################################################################################################################
## Make a dataframe with all of the appropriate information
##   The first column must contain the filenames of the configs that will be generated
##   The second column must contain the jobnames that will be displayed as each file is running
##   All other columns must contain the information to be inserted into the config files, with the column headers matching the
##     appropriate line number in the config file.

#' USDOS: Generate config files
#' 
#' A streamlined method for producing large numbers of USDOS config files and changing parameters. For additional details, refer to the USDOS User Manual. 
#' 
#' @param run.control The initial parameter to assign the control to be implemented. This value changes the default values for specific additional parameters. The default value is 'noControl', indicating no controls will be implemented. Other options include "MB" (movement ban), 'MB_IPcull' (movement ban + infected premises (IP) culling), 'MB_IPDCcull (movement ban  + culling of IPs and dangerous contacts (DCs))', 'MB_cullVax' (movement ban, IP culling, and vaccination either of DCs or in a ring), 'IPcull' (IP culling), and "other". For "other" runs, the user must specify all control-related inputs manually. 
#' @param parameter.sample Used only for sensitivity runs. Options are "spread" or "control". "spread" is used to provide variation in the tranmissibility constant (transm.const), K2, K3, latency, and infectious parameters, using values provided in a supplemental file called "LHS_Parameters.csv" in the inputfiles folder. "control" is used to provide variation in cull.rate, vax.rate, effective.vax, index.rep.time, rep.time, DC.rep.time, DC.sus, and DC.exp parameters, using values provided in a supplemental file called "LHS_Control_Parameters_040119.csv" in the inputfiles folder. The file must contain a number of rows equal to the total number of simulations being generated. Default is NA (a non-sensitivity run).
#' @param sens.number Used only for sensitivity runs.The row number in the sensitivity parameter file to use.
#' @param batch.name (config line 1) Batch name, used as prefix for output files. Only letters, numbers and underscore. If not specified, a name will be generated based on run options. 
#' @param summary (2) Generate a file with summary information about each seeded outbreak (0 = no, 1 = yes (default))
#' @param detail (3) Generate a detail file withinformation on each exposure event (0 = no, 1 = yes (default))
#' @param print.grid (4) Print grid cells on/off (0 = no (default), 1 = yes )
#' @param flaps (11) The subname of the FLAPS version being used. For example, 'flaps12_min' will use the files with a naming pattern 'flaps12_min_00XX.txt'
#' @param species (12) List of species for which counts are provided in premises file, comma-separated. Default is "beef,dairy".
#' @param timesteps (13) The number of timesteps to run for each seeded outbreak. The timestep is days for foot and mouth disease (FMD).
#' @param cutoff (14) The maximum number of infected premises a seeded outbreak may generate. If this is reached, the run will be terminated. '*' indicates no limit (default).
#' @param verbose (15) Output extra information to console as model runs: 0 = off, 1 = basic steps (default), 2 = detailed information for debugging
#' @param xy (17) Reverse x/y option: if input file is entered as lat/long (y/x) set to 1, and if long/lat (x/y) set to 0 (default)
#' @param fips.info (18) Name of file containing fips name, state name, area (m2), x, y. Tab separated. Default is "inputfiles/FIPS_updated_new.txt"
#' @param start.day (19) Day of the year to start the simulation (Jan 1 = day 1). Must be 0 - 365. 0 = random day [1-365] (default). 
#' @param disease (20) The type of disease to simulate: 0 = foot and mouth disease. 
#' @param source (21) The source FIPS to be seeded. Either a filename or 'allFips' (default). 'allFips' will seed in all counties containing premises. Filename is of the tab-delimited file containing the FIPS codes or premisesIDs from which to seed infection, with one line per simulation. To seed from multiple premises at once, provide comma-separated premisesIDs in the file. 
#' @param filetype (22) The type of information provided in the file entered for "source". 'fips' (default) will choose 1 premises at random within the FIPS codes listed in the source parameter or in all FIPS, 'singlePremises' will use the premisesID provided by the source parameter, 'multiplePremises' will seed all comma separated premisesIDs in the source parameter.
#' @param market.within (23) Controls market behavior with regard to within-herd spread. For FMD, this is a real number probability (range [0,1], default = 0.5) that a shipment from an infected market transmits infection. 
#' @param susc.exp.sp1 (24) Susceptibility exponent (power, q) for the first species listed on line 12 of the config file. Default species is beef, with power = 1. Only affects FMD runs. 
#' @param susc.exp.sp2 (24) Susceptibility exponent (power, q) for the second species listed on line 12 of the config file. Default species is dairy, with power = 1. Only affects FMD runs.  
#' @param transm.exp.sp1 (25) Transmissibility (infectiousness) exponent (power, p) for the first species listed on line 12 of the config file. Default species is beef, with power = 1. Only affects FMD runs.  
#' @param transm.exp.sp2 (25) Transmissibility (infectiousness) exponent (power, p) for the second species listed on line 12 of the config file. Default species is dairy, with power = 1. Only affects FMD runs.  
#' @param susc.const.sp1 (26) Susceptibility constant for the first species listed on line 12 of the config file. Default species is beef, with constant = 1. Only affects FMD runs. 
#' @param susc.const.sp2 (26) Susceptibility constant for the second species listed on line 12 of the config file. Default species is dairy, with constant = 1. Only affects FMD runs. 
#' @param transm.const.sp1 (27) Transmissibility (infectiousness) constant for the first species listed on line 12 of the config file. Default species is beef, with constant = 10.252. Only affects FMD runs.  
#' @param transm.const.sp2 (27) Transmissibility (infectiousness) constant for the second species listed on line 12 of the config file. Default species is dairy,  with constant = 10.252. Only affects FMD runs.  
#' @param kernel.type (28) Kernel type for local (diffusion) spread: 0: k1/(1 + (distance/k2)^k3), 1: data file provided with data.kernel.file , 2: k1/(1+d/k2)^k3. Only affects FMD runs. 
#' @param k1 (29) Kernel parameter 1. Default is 1.46e-08. Only affects FMD runs. 
#' @param k2 (29) Kernel parameter 2. Default is 1686.155. Only affects FMD runs. 
#' @param k3 (29) Kernel parameter 3. Default is 2.267. Only affects FMD runs. 
#' @param data.kernel.file (30) Name of file containing data-based local spread probabilities by distance, ex inputfiles/UKDataKernel.txt. Only affects FMD runs. 
#' @param latency (31) The mean days from premises exposure to infectiousness. Default is 5 days. Only affects FMD runs. 
#' @param latency.sd (31) The standard deviation days from premises exposure to infectiousness. Default is 0 days. Only affects FMD runs. 
#' @param infectious (32) The mean days from premises infectiousness to immunity. Default is 20 days. Only affects FMD runs. 
#' @param infectious.sd (32) The standard deviation days from premises infectiousness to immunity. Default is 0 days. Only affects FMD runs. 
#' @param partial (33) The partial transition flag (0 = off, 1 = on (default)). Only affects FMD runs. 
#' @param partial.param (34) Partial transition parameters (of which there are 5). Either enter all 5 (comma separated) or "*" if not using partial transition. Only affects FMD runs. 
#' @param grid.file (36) Filename containing grid cells (will override other grid-related options)
#' @param grid.side (37) Length of cell side for uniform cells (will override grid density options grid.max.farms and grid.min.side)
#' @param grid.max.farms (38) Max farms per grid cell (default 500)
#' @param grid.min.side (38) Grid cell size minimum side length in meters (default 100,000)
#' @param shipment.gen (41) Method to generate county-to-county shipments with USAMM, comma-separated (0: shipments off; 1: USAMMv1; 3: USAMMv2 kernel 1; 4: USAMMv2 kernel 2; 5: USAMMv2 kernel 3 (default))
#' @param shipment.scale (42) Shipment scaling factor (default 1)
#' @param usamm.posteriors (44) USAMM posterior files. Must match method selected in shipment.gen. Comma-separated, one file for each species. Default is inputfiles/beef_k3_cov.post, inputfiles/dairy_k3_cov.post.
#' @param usamm.order (45) The order in which the temporal switching of USAMM parameters should happen. Comma separated, time periods exactly as the temporal component of the USAMM parameter names in the usamm.posterior file(s). Default is Q1,Q2,Q3,Q4.
#' @param usamm.day (46) Day of the year (Jan 1 = day 1) where each time period in usamm.order begins. Comma separated integers. Assume no leap years. Default is 1,91,182,274.
#' @param usamm.origin.cov (47) Origin covariates for all counties. One file for each species, comma separated. Must have header: FIPS, name1, name2... Names must match name component of covariate parameters in USAMM parameter file. Default is inputfiles/county_covariates_scaled_mean_zero.txt, inputfiles/county_covariates_scaled_mean_zero.txt.
#' @param usamm.dest.cov (48) Destination covariates for all counties. One file for each species, comma separated. Must have header: FIPS, name1, name2... Names must match name component of covariate parameters in USAMM parameter file. Default is inputfiles/county_covariates_scaled_mean_zero.txt, inputfiles/county_covariates_scaled_mean_zero.txt.
#' @param exp.ship (50) Exposed shipments: shipping from a farm with status exposed will cause the receiver to become exposed as well (0 = off, 1 = on (default)). Affects FMD runs only.
#' @param ctrl.type (51) The control type(s) implemented. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param ctrl.constraint.type (52) Control constraint functions. Options are noLimit, dailyLimit, stateSum (cull only), and nationalLimit (vax only). Set by default with run.control.
#' @param ctrl.constraint (53) Control constraint parameters, comma-separated parameters, SEMICOLON-separated by type. Movement bans have a single value (which is 0 by default), whereas other controls require two values: a mean and standard deviation. Set by default with run.control. Required when run.control = "other".
#' @param cull.rate (53) The mean number of animals per timestep that can be culled on a premises. Set by default with run.control. Required when run.control = "other".
#' @param cull.rate.sd (53) The standard deviation number of animals per timestep that can be culled on a premises. Set by default with run.control. Required when run.control = "other".
#' @param vax.rate (53) The mean number of animals per timestep that can be vaccinated on a premises. Set by default with run.control. Required when run.control = "other".
#' @param vax.rate.sd (53) The standard deviation number of animals per timestep that can be vaccinated on a premises. Set by default with run.control. Required when run.control = "other".
#' @param ctrl.scale (54) Spatial scale at which control is applied, comma-separated (fixed options: premises, county (shipBan only), state (shipBan only)). Set by default with run.control. Required when run.control = "other".
#' @param ctrl.constraint.files (55) List of files containing control constraints. Only used when run.control= "other". 
#' @param landfill.files (55) Set to TRUE if there is a unique landfill file for each FLAPS file (default). If landfill files are not named "landfills_flaps12_{flaps number}.txt" (e.g. 'landfills_flaps12_0001.txt'), specify naming convention with landfill.file.name. Otherwise, set to FALSE and the generic 'inputfiles/landfills_formatted.txt' file will be used. Set by default with run.control. 
#' @param landfill.file.name (55) Specify prefix for the landfill filenames (i.e. the part before "_{FLAPS number}.txt"). By default, "landfills_flaps12". 
#' @param vaccine.file (55) Set to the filename if there is a vaccination constraint file. The default for when vaccination is on will be a vaccination file named "vaccineBankUpdated.txt" in the inputfiles folder. Otherwise, default is "NA" for no vaccination file.
#' @param ctrl.constraint.filetypes (56) List of filetypes for control constraints listed in ctrl.constraint.files. Only used when run.control= "other"
#' @param effective.mean (57) Implemented to effective: mean number of timesteps, comma-separated for each control type. Only used for run.control ="other"
#' @param effective.mb (57) Implemented to effective (mean number of timesteps) for a movement ban. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param effective.cull (57) Implemented to effective (mean number of timesteps) for culling. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param effective.vax (57) Implemented to effective (mean number of timesteps) for vaccination. The timesteps between vaccination and vaccine effectiveness. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param effective.sd (58) Implemented to effective: standard deviation number of timesteps, comma-separated for each control type. Only used for run.control ="other"
#' @param effective.mb.sd (58) Implemented to effective (standard deviation number of timesteps) for a movement ban. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param effective.cull.sd (58) Implemented to effective (standard deviation number of timesteps) for culling. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param effective.vax.sd (58) Implemented to effective (standard deviation number of timesteps) for vaccination. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param inactive.mean (59) Effective to inactive: mean number of timesteps, comma-separated for each control type. Only used for run.control ="other".
#' @param inactive.mb (59) Effective to inactive (mean number of timesteps) for a movement ban. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param inactive.cull (59) Effective to inactive (mean number of timesteps) for culling. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param inactive.vax (59) Effective to inactive (mean number of timesteps) for vaccination. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param inactive.sd (60) Effective to inactive: standard deviation number of timesteps, comma-separated for each control type. Only used for run.control ="other".
#' @param inactive.mb.sd (60) Effective to inactive (standard deviation number of timesteps) for a movement ban, Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param inactive.cull.sd (60) Effective to inactive (standard deviation number of timesteps) for culling. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param inactive.vax.sd (60) Effective to inactive (standard deviation number of timesteps) for vaccination. Set by default with run.control, but can be specified using this option when run.control = "other". 
#' @param ctrl.eff (61) Control effectiveness (including compliance) of control types as proportion, comma-separated probability of preventing exposure given exposure and probability of transmission given infectiousness. SEMICOLON-separated by type. Only used for run.control="other"
#' @param mb.eff (61) The effectiveness of the movement ban put into place on shipments. The value is between 0 and 1. Set by default with run.control.
#' @param cull.eff (61) The effectiveness of culling. The value is between 0 and 1. Set by default with run.control. 
#' @param vax.eff (61) The effectiveness of the vaccination. The value is between 0 and 1. Set by default with run.control.
#' @param ctrl.triggers (63) Triggers for control implementation, comma-separated (fixed options: newRegionReportsOverX,newPremReportsOverX). '*' will turn off control. Set by default with run.control, required for run.control = "other".
#' @param ctrl.trigger.threshold (64) Control trigger thresholds, comma-separated (numeric). How many of "ctrl.triggers" must occur for control to be implemented? Set by default with run.control, required for run.control = "other".
#' @param ctrl.trigger.response (65) What measures are taken in response to meeting the control trigger threshold? Must exist in ctrl.type. Comma-separated. Set by default with run.control, required for run.control = "other".
#' @param ctrl.response.target (66) The response targets, comma-separated. -1 = DCs, 0 = triggers only, # = radius in units of x/y coordinates of premises. Only used if run.control="other"
#' @param vax.range (66) The range of the vaccination ring. The default value is '-1', indicating that the vaccination is applied to dangerous contacts, and not in a ring. A positive number indicates a ring vaccination, in meters. This is only in effect if a vaccination control type is selected.
#' @param ctrl.response.priority (67) Method for prioritizing premises for control. Currently a fixed option of "earliest" for all controls. 
#' @param index.rep.time (71) The mean timesteps from index premises exposure to reporting. Default is 15.
#' @param index.rep.time.sd (71) The standard deviation timesteps from index premises exposure to reporting.  Default is 0 for both diseases.
#' @param rep.time (72) The mean timesteps from non-dangerous contact premises exposure to reporting. Default is 8.
#' @param rep.time.sd (72) The standard deviation timesteps from non-dangerous contact premises exposure to reporting. Default is 0 for both diseases.  
#' @param DC.rep.time (73) The mean timesteps from dangerous contact premises exposure to reporting. Default is 2.
#' @param DC.rep.time.sd (73) The standard deviation timesteps from dangerous contact premises exposure to reporting. Default is 0 for both diseases.
#' @param DC.sus (74) The scale for tracing DCs of a susceptible premises. The default value is 4 to make the mean DCs per reported premises between 1.5-2 when combined with the exposed value.
#' @param DC.exp (74) The scale for tracing DCs of a exposed premises. The default value is 5 to make the mean DCs per reported premises between 1.5-2 when combined with the susceptible value.
#' @param template.config.location The location of the template for the config file. Default is "templates/USDOS_CONFIG_TEMPLATE_General.txt".
#' @param sh.template The destination of the SH_FILE_TEMPLATE.txt file. Default is 'templates/SH_FILE_TEMPLATE.txt'.
#' @param jobfile.name The name of the jobfile to output. Default is 'USDOS_Sim.job'. 
#' @param destination.folder The location to put all generated files. Folder must exist. Default is the current working directory.
#' @param replacement.df The dataframe containing all information to use to modify the default config file. This is created with the time, date, and replacement_df as its title.
#' @param slurm If running on a SLURM system, set option to 1 to generate BATCH files and a jobfile that lists BATCH file names. Setting to 0 (default) indicates files will not be used on a SLURM system. BATCH files will not be genreated and the jobfile will reference config file names. 
#' @param reps The number of times each FLAPS will be repeated. The default value is 10, for a total of 100 runs. This is the minimum number of runs that should be interpreted to adequately capture uncertainty. 

#' @return This code will generate configuration (config) files and other necessary files to run USDOS (job and BATCH files), as well as a csv of the replacement data frame with the options used to fill in the configs.


#' @export

##' @author Send bug reports, suggestions, corrections, or comments to webblaboratory(at)gmail.com
##' 
######################################################
## Start Function ##

createConfigs <- function(run.control = "noControl",  disease = 0, ...){ 
  
  options(scipen=3)  
  
  ### Set default values for config files: 
  
  # Inputs
  template.config.location <- "templates/USDOS_CONFIG_TEMPLATE_General.txt"
  sh.template <- "templates/SH_FILE_TEMPLATE.txt"
  jobfile.name <- "USDOS_Sim.job"
  batch.name <- ""  
  slurm <- 0
  
  # Outputs
  destination.folder <- getwd()
  replacement.df <- data.frame()
  
  # General
  summary <- 1
  detail <- 1
  print.grid <- 0
  
  flaps <- "FLAPS12_Quarterly_USDOS_format"
  species <- "beef,dairy"
  
  verbose <- 1
  reps <- 10
  cutoff <- "*"  
  
  fips.info <- "inputfiles/FIPS_updated_new.txt"  
  xy <- 0
  start.day <- 0
  
  source <- "allFips"
  filetype <- "fips"
  
  kernel.type <- 0
  data.kernel.file <- "*"
  
  landfill.files <- TRUE
  landfill.file.name <- "landfills_flaps12"
  
  # grid
  grid.file <- "*"
  grid.side <- "*"
  grid.max.farms <-500
  grid.min.side <- 100000
  
  # Shipment
  shipment.gen <- 5
  shipment.scale <- 1
  usamm.posteriors <- "inputfiles/beef_k3_cov.post, inputfiles/dairy_k3_cov.post"
  usamm.order <- "Q1,Q2,Q3,Q4"
  usamm.day <- "1,91,182,274"
  usamm.origin.cov <- "inputfiles/county_covariates_scaled_mean_zero.txt,inputfiles/county_covariates_scaled_mean_zero.txt"
  usamm.dest.cov <-   "inputfiles/county_covariates_scaled_mean_zero.txt,inputfiles/county_covariates_scaled_mean_zero.txt"
  
  exp.ship <-1 # only affects FMD
  
  # Sensitivity 
  parameter.sample <- "NA"
  sens.number <- NA
  
  # Disease settings. 
  DC.sus <- 4
  DC.exp <- 5  
  index.rep.time.sd <- 0
  rep.time.sd <- 0
  DC.rep.time.sd <- 0
  cull.rate.sd <-0
  

  k1 <- 1.46e-08
  k2 <- 1686.155
  k3 <- 2.267
  latency <- 5
  latency.sd <- 0
  infectious <- 20 
  infectious.sd <- 0
  susc.exp.sp1 <- 1
  susc.exp.sp2 <- 1
  transm.exp.sp1 <- 1
  transm.exp.sp2 <- 1
  susc.const.sp1 <- 1
  susc.const.sp2 <- 1
  transm.const.sp1 <- 10.252
  transm.const.sp2 <- 10.252
  
  # Partial Transition
  partial <- 1
  partial.param <- "0.05,0.006,0.44,4,6.30852"
  

  
  # Settings that need to be disease-specific
  if(disease == 0){ # FMD
    timesteps <- 365
    market.within <- 0.5
    index.rep.time <- 15
    rep.time <- 8
    DC.rep.time <- 2
    cull.rate <- 240 } 
  
  
  ### Default control settings for a base (noControl) run
  ctrl.type <- "*"      
  mb.eff <- 0
  cull.eff <- 1
  ctrl.constraint.type <-  "noLimit"
  ctrl.constraint <- 0 
  ctrl.scale <- "state"
  ctrl.constraint.files <- NULL
  ctrl.constraint.filetypes <- NULL
  
  ctrl.triggers <- "*" 
  ctrl.trigger.threshold <- 0
  ctrl.trigger.response <-"*"
  ctrl.response.target <- 0
  ctrl.response.priority <- "earliest"
  
  effective.mean <-NULL
  effective.sd <- NULL
  effective.mb <- 0
  effective.mb.sd <- 0
  effective.cull <- 0
  effective.cull.sd <- 0
  effective.vax <- 11 
  effective.vax.sd <- 0
  
  inactive.mean <- NULL
  inactive.sd <- NULL
  inactive.mb <- 366
  inactive.mb.sd <- 0
  inactive.cull <- 366
  inactive.cull.sd <- 0
  inactive.vax <- 183
  inactive.vax.sd <- 0
  
  vax.rate <- NULL
  vax.rate.sd <- NULL
  vax.eff <- 0
  vax.range <- NULL
  vaccine.file <- "NA"
  
  
  ################################################################################################################
  
  ### Set the control type with if statements. ####
  
  # Movement ban only 
  if(run.control == "MB"){
    ctrl.type <- "shipBan"
    ctrl.constraint.type <-  "noLimit"
    ctrl.scale <- "state"
    mb.eff <- 0.75 
    ctrl.triggers <- "newRegionReportsOverX"
    ctrl.trigger.threshold <- 0
    ctrl.trigger.response <- "shipBan"
  }
  
  
  if(run.control == "MB_IPcull"){
    ctrl.type <- "shipBan,cull"
    ctrl.constraint.type <-  "noLimit,stateSum"
    ctrl.scale <- "state,premises"
    mb.eff <- 0.75 
    cull.eff <- 1
    ctrl.triggers <- "newRegionReportsOverX,newPremReportsOverX"
    ctrl.trigger.threshold <- "0,0"
    ctrl.trigger.response <- "shipBan,cull"
  }
  
  
  if(run.control == "MB_IPDCcull"){
    ctrl.type <- "shipBan,cull"
    ctrl.constraint.type <-  "noLimit,stateSum"
    ctrl.scale <- "state,premises"
    mb.eff <- 0.75 
    cull.eff <- 1
    ctrl.triggers <- "newRegionReportsOverX,newPremReportsOverX,newPremReportsOverX"
    ctrl.trigger.threshold <- "0,0,0"
    ctrl.trigger.response <- "shipBan,cull,cull"
  }
  
  if(run.control == "MB_cullVax"){
    ctrl.type <- "shipBan,cull,vax"
    ctrl.constraint.type <-  "noLimit,stateSum,nationalLimit"
    vax.rate <- 6804 
    vax.rate.sd <- 0
    ctrl.scale <- "state,premises,premises"
    mb.eff <- 0.75 
    cull.eff <- 1
    vax.eff <- 0.9
    vax.range <- -1 # default to DC vaccination
    ctrl.triggers <- "newRegionReportsOverX,newPremReportsOverX,newPremReportsOverX"
    ctrl.trigger.threshold <- "0,0,0"
    ctrl.trigger.response <- "shipBan,cull,vax"
    vaccine.file <-  "vaccineBankUpdated.txt"
  }
  
  if(run.control == "IPcull"){
    
    ctrl.type <- "cull"
    ctrl.constraint.type <-  "stateSum"
    ctrl.scale <- "premises"
    cull.eff <- 1
    ctrl.triggers <- "newRegionReportsOverX"
    ctrl.trigger.threshold <- "0"
    ctrl.trigger.response <- "cull"
  }
  
  
  
  
  ################################
  ## Input checks ####
  
  ## List of all allowed function inputs 
  additional.inputs = list(...)
  parameter.options = c("template.config.location", "sh.template", "jobfile.name", "summary","detail","print.grid","flaps", 
                        "species","xy","start.day","landfill.files", "landfill.file.name","vaccination.file",
                        "fips.info", "reps", "cutoff", "timesteps", "source", "filetype", "parameter.sample", 
                        "k1", "k2", "k3","susc.exp.sp1","susc.exp.sp2","transm.exp.sp1", "susc.const.sp1","susc.const.sp2","transm.const.sp1","transm.const.sp2", 
                        "transm.exp.sp2","kernel.type","data.kernel.file","latency",'latency.sd', "infectious", "infectious.sd",
                        "grid.file","grid.side","grid.max.farms","grid.min.size",
                        "shipment.gen","shipment.scale","usamm.posteriors","usamm.order","usamm.day","usamm.origin.cov","usamm.dest.cov",
                        "mb.eff", "cull.rate", "cull.rate.sd","vax.rate","vax.rate.sd", "ctrl_type", 
                        "ctrl.constraint.files","ctrl.constraint.filetypes", "cull.eff", "ctrl.eff",
                        "effective.mean", "effective.sd" , "effective.mb" , "effective.mb.sd" , "effective.cull" , "effective.cull.sd" , "effective.vax" ,
                        "effective.vax.sd","inactive.mean","inactive.sd","inactive.mb" , "inactive.mb.sd" , "inactive.cull" , "inactive.cull.sd" ,
                        "inactive.vax" ,"inactive.vax.sd","vax.eff",  "vax.range", "partial.param", "partial", "index.rep.time", "rep.time", 
                        "DC.rep.time","index.rep.time.sd", "rep.time.sd","DC.rep.time.sd", "DC.sus", "DC.exp",
                        "destination.folder", "replacement.df", "sens.number","exp.ship","verbose",
                        "batch.name","market.within", "ctrl.type",
                        "ctrl.constraint.type","ctrl.constraint","ctrl.scale","ctrl.triggers","ctrl.trigger.threshold",
                        "ctrl.trigger.response","ctrl.response.target", "ctrl.response.priority","slurm")
  
  ## Return error if additional inputs aren't potential parameters
  unknown.inputs = additional.inputs[!names(additional.inputs) %in% parameter.options]
  if(length(unknown.inputs) > 0){
    error.message = paste0("----------------------------------------------------\n",
                           "These arguments aren't recognized: \n",
                           paste(names(unknown.inputs), collapse = "\n"), "\n",
                           "Files not generated. \n",
                           "----------------------------------------------------\n")
    return(cat(error.message))
  }
  ## Use only the additional arguments that match the possible parameter names
  additional.inputs = additional.inputs[names(additional.inputs) %in% parameter.options]
  
  ## Override values based on additional arguments passed to the function. 
  for(parameter in names(additional.inputs)){
    assign(parameter, additional.inputs[parameter][[1]]) 
  }
  
  ## For sensitivity runs, create parameter values based on the latin hypercube sampling ####
  if(parameter.sample == "spread"){
    lh <- read.csv("inputfiles/LHS_Parameters.csv", header = TRUE)
    transm.const.sp1  <- lh$Tc[sens.number]
    k1 <- lh$k1[sens.number]
    k2 <- lh$k2[sens.number]
    k3 <- lh$k3[sens.number]
    latency <- lh$lat[sens.number]
    infectious <- lh$inf[sens.number]
  }
  
  if(parameter.sample == "control"){
    lh <- read.csv("inputfiles/LHS_Control_Parameters_040119.csv", header = TRUE)
    cull.rate <- lh$cull_rate[sens.number]
    vax.rate <- lh$vax_rate[sens.number]
    effective.vax <- lh$vax_delay[sens.number]
    index.rep.time <- lh$index_rep_time[sens.number]
    rep.time <- lh$rep_time[sens.number]
    DC.rep.time <- lh$DC_rep_time[sens.number]
    DC.sus <- lh$DC_sus[sens.number]
    DC.exp <- lh$DC_exp[sens.number]
    mb.eff <- lh$shipban[sens.number]
    if(run.control == "IPDCcull" | run.control == "IPcull"|run.control == "MB_IPDCcull" | run.control == "MB_IPcull"){vax.rate = NULL}
    if(run.control == "IPDCcull" | run.control == "IPcull"|run.control == "MB_IPDCcull" | run.control == "MB_IPcull"){effective.vax = NULL}
  } 
  
  ##########################################
  ## Creating FLAPS names and Landfill files ####
  ## Create list of FLAPS files from the inputted FLAPS name ("flaps" argument)
  flaps.names <- paste0(flaps,"_",  sprintf("%04d", 1:10),".txt", sep = "")
  all.flaps.names <- rep(flaps.names, each = reps)
  
  ## Create list of additional control files from the inputted control argument  
  if(landfill.files == TRUE){
    landfill.file.list <- paste0(landfill.file.name,"_", sprintf("%04d", 1:10), ".txt", sep = "")
    landfill.file.names <- rep(landfill.file.list, each = reps)} else {
      landfill.file.names <- rep('landfills_formatted.txt', each = 10 * reps)}
  
  ##########################################
  ## Create file naming conventions #### 
  
  # vaccine ring radius (also includes check for allowable values of vax.range)
  vax_rng <-if(is.null(vax.range) == T){-1}else if(vax.range == -1){-1} else if(vax.range>=0){vax.range} else stop ("Vaccination Range must be positive or '-1' if not using ring vaccination")
  
  vax <- if(vax_rng>=0){paste0(vax_rng/1000,"km_")}else{""} 
  
  # Movement ban effectiveness
  MvmtBan <- if(mb.eff!=0){paste0("MvmtBan_",100*mb.eff,"_")}else{""}
  
  # Shipments off
  ShipOff <-if(shipment.gen == 0) {"ShipmentsOff_"} else{""}
  
  # Errors for parameters that are out of bounds and warnings for unique situations
  if(mb.eff < 0 | mb.eff > 1) stop ('Movement Ban Effectiveness (mb.eff) must be between 0 and 1')
  if(vax.eff < 0 | vax.eff > 1) stop ('Vaccine Effectiveness (vax.eff) must be between 0 and 1')
  
  # Generate batch.name if one isn't provided. Add pieces to output filenames here (ex for FMD, whether partial transition is on and the ifnectous period).
  # This batch.name is what's used in post-processing to determine unique run types
  batch.name = ifelse(batch.name =="",paste0("FMD_",ifelse(disease==0,paste0(ifelse(partial==1,"PTon","PToff"),"_Infectious",infectious,"days")),"_",
                                             run.control,"_",ifelse(vax == "","",vax),ifelse(MvmtBan == "","",MvmtBan),ifelse(ShipOff == "","",ShipOff)),
                      batch.name)
  
  ## Create the replacement_df data frame that contains all of the information created by this function.  ####
  config.fname <- paste0("config_", batch.name, "_",
                         gsub(".txt", "", all.flaps.names), "_", sprintf("%02d", rep(1:10, 10)), 
                         format(Sys.time(), '_%Y%m%d'), ifelse(is.na(sens.number), "", paste0("_",sens.number)), ".txt")
  
  jobname <- paste0("MI_",batch.name,"_f", sprintf("%02d", rep(1:10, each=10)), "_", 
                    sprintf("%02d", rep(1:10),"_", sens.number), format(Sys.time(), '_%Y%m%d'))
  
  replacement.df <- data.frame(config.fname, jobname, 
                               summary = summary,
                               detail = detail,
                               print.grid = print.grid,
                               FLAPS.locations = paste0("FLAPS/", all.flaps.names), 
                               species = species,
                               timesteps = paste0(timesteps),
                               cutoff = paste0(cutoff),
                               verbose=verbose,
                               xy = xy,
                               fips.info = paste0(fips.info),
                               start.day = start.day,
                               disease = disease,
                               source = paste0(source), 
                               filetype = paste0(filetype),
                               market.within=market.within,
                               susceptibility.exponents=paste0(susc.exp.sp1,",",susc.exp.sp2),
                               infectiousness.exponents=paste0(transm.exp.sp1,",",transm.exp.sp2),
                               susceptibility.constants=paste0(susc.const.sp1,",",susc.const.sp2),
                               infectiousness.constants=paste0(transm.const.sp1,",",transm.const.sp2),
                               kernel.type = kernel.type,
                               kernels = paste0(k1,",", k2, ",", k3), 
                               data.kernel.file = data.kernel.file,
                               latency.dist= paste0(latency, ",",latency.sd), 
                               infectious.dist = paste0(infectious, ",",infectious.sd),
                               partial = partial, 
                               partial.param = partial.param,
                               grid.file = grid.file,
                               grid.side = grid.side,
                               grid.params = paste0(grid.max.farms,",",grid.min.side),
                               shipment.gen = paste0(shipment.gen),
                               shipment.scale = paste0(shipment.scale),
                               usamm.posteriors = paste0(usamm.posteriors),
                               usamm.order = paste0(usamm.order),
                               usamm.day = paste0(usamm.day),
                               usamm.origin.cov = paste0(usamm.origin.cov),
                               usamm.dest.cov = paste0(usamm.dest.cov),
                               exp.ship=exp.ship,
                               ctrl.type = ctrl.type,
                               ctrl.constraint.type = ctrl.constraint.type,
                               ctrl.constraint = if(run.control == "other") {ctrl.constraint <- ctrl.constraint} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {0} else
                                     if (run.control == "MB_IPcull" |run.control == "MB_IPDCcull"){paste0("0;",cull.rate,",",cull.rate.sd)} else
                                       if (run.control == "MB_cullVax" ){paste0("0;",cull.rate,",",cull.rate.sd,";",vax.rate,",",vax.rate.sd)} else
                                         if (run.control == "IPcull") {paste0(cull.rate,",",cull.rate.sd)},
                               ctrl.scale = ctrl.scale,
                               ctrl.constraint.files = if(run.control == "other") {ctrl.constraint.files} else 
                                 if (run.control == "noControl") {"NA"} else
                                   if (run.control == "MB") {"NA"} else
                                     if (run.control == "MB_IPcull" |run.control == "MB_IPDCcull"){paste0("NA;inputfiles/",landfill.file.names)} else
                                       if (run.control == "MB_cullVax"){paste0("NA;inputfiles/",landfill.file.names,";inputfiles/",vaccine.file)} else
                                         if (run.control == "IPcull") {paste0("inputfiles/", landfill.file.names)},
                               ctrl.constraint.filetypes = if(run.control == "other") {ctrl.constraint.filetypes} else 
                                 if (run.control == "noControl") {"NA"} else
                                   if (run.control == "MB") {"NA"} else
                                     if (run.control == "MB_IPcull" |run.control == "MB_IPDCcull" ){"NA;resourceLocs"} else
                                       if (run.control == "MB_cullVax"){"NA;resourceLocs;resourceBoosts"} else
                                         if (run.control == "IPcull") {"resourceLocs"},
                               effective.mean = if(run.control == "other") {effective.mean} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {effective.mb} else
                                     if (run.control == "MB_IPcull" ) {paste0(effective.mb,",",effective.cull)} else
                                       if (run.control == "MB_IPDCcull") {paste0(effective.mb,",",effective.cull)} else
                                         if (run.control == "MB_cullVax") {paste0(effective.mb,",",effective.cull,",",effective.vax)} else
                                           if (run.control == "IPcull") {effective.cull},
                               effective.sd = if(run.control == "other") {effective.mean.sd} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {effective.mb.sd} else
                                     if (run.control == "MB_IPcull" ) {paste0(effective.mb.sd,",",effective.cull.sd)} else
                                       if (run.control == "MB_IPDCcull") {paste0(effective.mb.sd,",",effective.cull.sd)} else
                                         if (run.control == "MB_cullVax") {paste0(effective.mb.sd,",",effective.cull.sd,",",effective.vax.sd)} else
                                           if (run.control == "IPcull") {effective.cull.sd},
                               inactive.mean = if(run.control == "other") {inactive.mean} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {inactive.mb} else
                                     if (run.control == "MB_IPcull" ) {paste0(inactive.mb,",",inactive.cull)} else
                                       if (run.control == "MB_IPDCcull") {paste0(inactive.mb,",",inactive.cull)} else
                                         if (run.control == "MB_cullVax") {paste0(inactive.mb,",",inactive.cull,",",inactive.vax)} else
                                           if (run.control == "IPcull") {inactive.cull},
                               inactive.sd = if(run.control == "other") {inactive.mean.sd} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {inactive.mb.sd} else
                                     if (run.control == "MB_IPcull" ) {paste0(inactive.mb.sd,",",inactive.cull.sd)} else
                                       if (run.control == "MB_IPDCcull") {paste0(inactive.mb.sd,",",inactive.cull.sd)} else
                                         if (run.control == "MB_cullVax") {paste0(inactive.mb.sd,",",inactive.cull.sd,",",inactive.vax.sd)} else
                                           if (run.control == "IPcull") {inactive.cull.sd},
                               ctrl.response.target = if(run.control == "other") {ctrl.response.target} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {0} else
                                     if (run.control == "MB_IPcull" ) {"0,0"} else
                                       if (run.control == "MB_IPDCcull") {"0,0,-1"} else
                                         if (run.control == "MB_cullVax") {paste0("0,0,",vax.range)} else
                                           if (run.control == "IPcull") {"0"},
                               ctrl.triggers = ctrl.triggers,
                               ctrl.eff = if(run.control == "other") { ctrl.eff} else 
                                 if (run.control == "noControl") {0} else
                                   if (run.control == "MB") {paste0(mb.eff, ",", mb.eff)} else
                                     if (run.control == "MB_IPcull" |run.control == "MB_IPDCcull"){paste0(mb.eff, ",", mb.eff,";",cull.eff, ",", cull.eff)} else
                                       if (run.control == "MB_cullVax" ){paste0(mb.eff, ",", mb.eff,";",cull.eff, ",", cull.eff,";", vax.eff, ",", vax.eff)} else
                                         if (run.control == "IPcull") {paste0(cull.eff, ",", cull.eff)},
                               ctrl.trigger.threshold = ctrl.trigger.threshold,
                               ctrl.trigger.response = ctrl.trigger.response,
                               ctrl.response.priority = if (run.control == "other") {ctrl.response.priority} else
                                 if (run.control == "noControl" | run.control == "MB" | run.control == "IPcull") {"earliest"} else
                                   if (run.control == "MB_IPcull") {paste0("earliest",",","earliest")} else
                                     if (run.control == "MB_IPDCcull" | run.control == "MB_cullVax") {paste0( "earliest",",","earliest",",","earliest")},
                               index.rep.time = paste0(index.rep.time,",",index.rep.time.sd ),
                               rep.time = paste0(rep.time,",",rep.time.sd),
                               DC.rep.time = paste0(DC.rep.time, ",",DC.rep.time.sd),
                               DC_scaling = paste0("sus,", DC.sus, ";exp,", DC.exp),
                               stringsAsFactors = FALSE)
  
  replacement.df$full.batch.name = gsub("config_", "", replacement.df$config.fname)
  replacement.df$full.batch.name = gsub(".txt", "", replacement.df$full.batch.name)
  
  columns <- length(replacement.df)
  for(i in 3:columns){
    if (colnames(replacement.df)[i] == "full.batch.name"){colnames(replacement.df)[i] = 1} else 
      if (colnames(replacement.df)[i] == "summary"){colnames(replacement.df)[i] = 2} else 
        if (colnames(replacement.df)[i] == "detail"){colnames(replacement.df)[i] = 3} else 
          if (colnames(replacement.df)[i] == "print.grid"){colnames(replacement.df)[i] = 4} else 
            if (colnames(replacement.df)[i] == "FLAPS.locations"){colnames(replacement.df)[i] = 11} else 
              if (colnames(replacement.df)[i] == "species"){colnames(replacement.df)[i] = 12} else 
                if (colnames(replacement.df)[i] == "timesteps"){colnames(replacement.df)[i] = 13} else 
                  if (colnames(replacement.df)[i] == "cutoff"){colnames(replacement.df)[i] = 14} else 
                    if (colnames(replacement.df)[i] == "verbose"){colnames(replacement.df)[i] = 15} else 
                      if (colnames(replacement.df)[i] == "xy"){colnames(replacement.df)[i] = 17} else
                        if (colnames(replacement.df)[i] == "fips.info"){colnames(replacement.df)[i] = 18} else
                          if (colnames(replacement.df)[i] == "start.day"){colnames(replacement.df)[i] = 19} else
                            if (colnames(replacement.df)[i] == "disease"){colnames(replacement.df)[i] = 20} else
                              if (colnames(replacement.df)[i] == "source"){colnames(replacement.df)[i] = 21} else 
                                if (colnames(replacement.df)[i] == "filetype"){colnames(replacement.df)[i] = 22} else 
                                  if (colnames(replacement.df)[i] == "market.within"){colnames(replacement.df)[i] = 23} else
                                    if (colnames(replacement.df)[i] == "susceptibility.exponents"){colnames(replacement.df)[i] = 24} else
                                      if (colnames(replacement.df)[i] == "infectiousness.exponents"){colnames(replacement.df)[i] = 25} else
                                        if (colnames(replacement.df)[i] == "susceptibility.constants"){colnames(replacement.df)[i] = 26} else
                                          if (colnames(replacement.df)[i] == "infectiousness.constants"){colnames(replacement.df)[i] = 27} else
                                            if (colnames(replacement.df)[i] == "kernel.type"){colnames(replacement.df)[i] = 28} else
                                              if (colnames(replacement.df)[i] == "kernels"){colnames(replacement.df)[i] = 29} else
                                                if (colnames(replacement.df)[i] == "data.kernel.file"){colnames(replacement.df)[i] = 30} else
                                                  if (colnames(replacement.df)[i] == "latency.dist"){colnames(replacement.df)[i] = 31} else
                                                    if (colnames(replacement.df)[i] == "infectious.dist"){colnames(replacement.df)[i] = 32} else
                                                      if (colnames(replacement.df)[i] == "partial"){colnames(replacement.df)[i] = 33} else
                                                        if (colnames(replacement.df)[i] == "partial.param"){colnames(replacement.df)[i] = 34} else
                                                          if (colnames(replacement.df)[i] == "grid.file"){colnames(replacement.df)[i] = 36} else
                                                            if (colnames(replacement.df)[i] == "grid.side"){colnames(replacement.df)[i] = 37} else
                                                              if (colnames(replacement.df)[i] == "grid.params"){colnames(replacement.df)[i] = 38} else
                                                                if (colnames(replacement.df)[i] == "shipment.gen"){colnames(replacement.df)[i] = 41} else
                                                                  if (colnames(replacement.df)[i] == "shipment.scale"){colnames(replacement.df)[i] = 42} else
                                                                    if (colnames(replacement.df)[i] == "usamm.posteriors"){colnames(replacement.df)[i] = 44} else
                                                                      if (colnames(replacement.df)[i] == "usamm.order"){colnames(replacement.df)[i] = 45} else
                                                                        if (colnames(replacement.df)[i] == "usamm.day"){colnames(replacement.df)[i] = 46} else
                                                                          if (colnames(replacement.df)[i] == "usamm.origin.cov"){colnames(replacement.df)[i] = 47} else
                                                                            if (colnames(replacement.df)[i] == "usamm.dest.cov"){colnames(replacement.df)[i] = 48} else
                                                                              if (colnames(replacement.df)[i] == "exp.ship"){colnames(replacement.df)[i] = 50} else
                                                                                if (colnames(replacement.df)[i] == "ctrl.type"){colnames(replacement.df)[i] = 51} else
                                                                                  if (colnames(replacement.df)[i] == "ctrl.constraint.type"){colnames(replacement.df)[i] = 52} else
                                                                                    if (colnames(replacement.df)[i] == "ctrl.constraint"){colnames(replacement.df)[i] = 53} else
                                                                                      if (colnames(replacement.df)[i] == "ctrl.scale"){colnames(replacement.df)[i] = 54} else
                                                                                        if (colnames(replacement.df)[i] == "ctrl.constraint.files"){colnames(replacement.df)[i] = 55} else 
                                                                                          if (colnames(replacement.df)[i] == "ctrl.constraint.filetypes"){colnames(replacement.df)[i] = 56} else
                                                                                            if (colnames(replacement.df)[i] == "effective.mean"){colnames(replacement.df)[i] = 57} else
                                                                                              if (colnames(replacement.df)[i] == "effective.sd"){colnames(replacement.df)[i] = 58} else
                                                                                                if (colnames(replacement.df)[i] == "inactive.mean"){colnames(replacement.df)[i] = 59} else
                                                                                                  if (colnames(replacement.df)[i] == "inactive.sd"){colnames(replacement.df)[i] = 60} else
                                                                                                    if (colnames(replacement.df)[i] == "ctrl.eff"){colnames(replacement.df)[i] = 61} else
                                                                                                      if (colnames(replacement.df)[i] == "ctrl.triggers"){colnames(replacement.df)[i] = 63} else
                                                                                                        if (colnames(replacement.df)[i] == "ctrl.trigger.threshold"){colnames(replacement.df)[i] = 64} else
                                                                                                          if (colnames(replacement.df)[i] == "ctrl.trigger.response"){colnames(replacement.df)[i] = 65} else
                                                                                                            if (colnames(replacement.df)[i] == "ctrl.response.target"){colnames(replacement.df)[i] = 66} else
                                                                                                              if (colnames(replacement.df)[i] == "ctrl.response.priority"){colnames(replacement.df)[i] = 67} else
                                                                                                                if (colnames(replacement.df)[i] == "index.rep.time"){colnames(replacement.df)[i] = 71} else
                                                                                                                  if (colnames(replacement.df)[i] == "rep.time"){colnames(replacement.df)[i] = 72} else
                                                                                                                    if (colnames(replacement.df)[i] == "DC.rep.time"){colnames(replacement.df)[i] = 73} else
                                                                                                                      if (colnames(replacement.df)[i] == "DC_scaling"){colnames(replacement.df)[i] = 74} else
                                                                                                                          warning(paste0(colnames(replacement.df)[i]," has not been added to this function. Change the line 
                                                                                                                                                   value in the template file then rerun this function without that variable."))
  }
  
  ## Write out the data frame with the associated parameters ####
  print(head(replacement.df))
  str1 <- paste0('replacementDF_',batch.name,'.csv') 
  write.csv(replacement.df, paste0(format(Sys.time(),'%Y%m%d_%H%M_'), sub('\\..*', '', str1), ".csv"), row.names = FALSE)
  
  ## Generate the Config, Batch, and Job files.
  destination.folder = paste0(destination.folder, "/")
  config = readLines(template.config.location)
  sh.template = readLines(sh.template)
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
    
    if (slurm == 1){
      #Generate a batch (.sh) file for each config:
      sh.file = gsub("JOBNAME", jobname, sh.template)
      sh.file = gsub("CONFIGNAME", fname, sh.file)
      fname.sh = gsub(".txt", ".sh", fname)
      fname.sh = gsub("config", "BATCH",fname.sh)
      cat(sh.file, file = paste0(destination.folder, fname.sh), sep = "\n")
      
      ### Add a line to the jobfile so the Batch file will run:
      cat(paste0("sbatch ", fname.sh, "\n"), file=paste0(destination.folder, jobfile.name), append=TRUE) } else
        # make jobfile that refers to configs instead of BATCH files
        cat(paste0(fname, "\n"), file=paste0(destination.folder, jobfile.name), append=TRUE)
    
  }
}





####################################################################################################################################
####################################################################################################################################

## Example commands ##

# # # FMD IPcull DCvax
# createConfigs(run.control = "MB_cullVax",destination.folder="./Generated_configs/TestBucket") 
# 
# # FMD IPcull ring vax
# createConfigs(run.control = "MB_cullVax",destination.folder="./Generated_configs/TestBucket",vax.range = 3000) 
# 
# # Shipments off, base run 
# createConfigs(shipment.gen = 0,destination.folder="./Generated_configs/TestBucket") 
# 
# # Only movement ban
# createConfigs(run.control = "MB",destination.folder="./Generated_configs/TestBucket") 

# 
# 
# # Sensitivity
# for(i in 1:100){
# 
# createConfigs(run.control = "MB_IPDCcull", 
#               source = "inputfiles/SensCounty_shortlist_0219.txt", 
#               parameter.sample = "control", 
#               sens.number = i,
#               jobfile.name = paste0("USDOS_Sens_", i, ".job"), 
#               destination.folder="./Generated_configs/TestBucket")
# }
# 
# 
# 

