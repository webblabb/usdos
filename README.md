USDOSv2.1
===

## Overview of the processes executed in "main.cpp"

1. Load and set parameters from configuration file @ File_manager::readConfig
2. Record parameters in runlog.txt
3. Load premises @ Grid_manager::readFarms
4. Create grid @ Grid_manager::initiateGrid
5. Set control parameters and load control resources @ Control_manager
6. Load/generate list of initial infections
7. Loop through list of initial infections, each is seed per simulation

     1. Initialize Status Manager, set initial infections @ Status_manager
     2. Initialize Shipment_manager
     3. Initialize Grid_checker
     4. Loop through timesteps or until epidemic dies or maximum threshold reached
          1. Advance timestep
          2. Auto-advance disease statuses @ Status_manager
          3. Auto-advance control statuses @ Control_manager
          4. Evaluate local spread (and dangerous contacts) @ Grid_checker
          5. Generate shipments from infectious premises @ Shipment_manager
          6. Filter shipments prevented by shipping bans @ Status_manager::filter_shipments
          7. Filter exposure prevented by premises level control @ Status_manager::eval_exposure
          8. Check for control actions triggered by simulation state and add premises/regions to waitlists @ Status_manager::add_waitlistMembers and @ Control_manager::check_controlRules
          9. Move premises/regions off waitlists and implement control as resources allow @ Status_manager::add_implemented
          10. Print detailed output for this timestep
     5. Print summary for this replicate

## Configuration file guide

The line number in the config file is denoted by \(#), along with a **brief descriptive variable name** and description including possible values. 

#### Output settings (lines 1-10)
\(1) **Batch name** A string containing only letters, numbers and underscore, used as prefix for output files. This is generated automatically based on selected options.  The date and time of run start are appended to the end of the file name.

\(2) **Summary output** Whether to generate the summary result file, which contains one line per replicate with number infected, duration, seed info, and run time. 

* **1** = Summary output on [default]. File will be outputted in working directory with name [Batch name]\_[date]\_[time]\_summary.txt  
* **0** = Summary output off 

\(3) **Detailed output** Whether to generate the detail result file, which contains one line per exposed premises with time, source of infection, and route of infection. 

* **1** = Detailed output on [default]. File will be outputted in working directory with name [Batch name]\_[date]\_[time]\_detail.txt
* **0** = Detailed output off 

\(4) **Print grid cells** Whether to generate a file containing the grid cells used. 

* **1** = Grid cells output on. File will be outputted in working directory with name [Batch name]_cells.txt
* **0** = Grid cells output off [default]
<!-- Get the wording right for what the grid cells are -->

\(5) Additional functionality under development. Do not change default value [0].

\(6) Additional functionality under development. Do not change default value [0].

\(7) Unused

\(8) Unused

\(9) Unused

\(10) Unused

#### General settings (lines 11-20)

\(11) **Premises file (FLAPS)** Path and name of file containing premises ID, FIPS (county ID), x, y, and population for each species- ie flaps12_min_0001.txt. By default this file is located in the FLAPS folder, and its path would thus be FLAPS/flaps12_min_0001.txt.  Make sure there are no numbers in scientific notation (ex. 1e+05). This may happen with IDs in R-formatted files.

\(12) **Species** List of species for which counts are provided in premises file, comma-separated. Default is "beef, dairy"

\(13) **Timesteps (days) to run** Positive integer of maximum number of days to run each simulation. If an outbreak ends before this time, the simulation will end at that point instead. Default is 365.

\(14) **Max infectious premises** Positive integer of threshold number of infectious premises allowed. If the total number of premises that have reached infectious status (not necessarily at the same time) surpasses this number, the simulation will end. The default value of \* means there is no limit on the number of infectious premises

\(15) **Verbose level** How much extra information to output to the console as model runs 

* **0** = Output only critical information to console 
* **1** = Output basic progress to console [default]
* **2** = Output debug-level information to console

\(16) Disabled 

\(17) **Reverse x/y**  This option allows you to specify if input files list locations as latitude(y) before longitude (x) (instead of the default, which is the reverse). 

* **0** = Location input file for premises (line 11) lists longitude (x) before latitude (y) [default]
* **1** = Location input file for premises (line 11) lists latitude(y) before longitude (x)

\(18) **FIPS info** Name and path of of file containing fips name, state name, area (in square meters), x coordinate, and y coordinate. Tab separated. Default is inputfiles/FIPS_20151805.txt

\(19) **Day of the year to seed outbreak** January 1 = day 1. Must be an integer 0 - 365. The default of 0 selects one random day [1-365].

\(20) Unused

#### Infection-related settings (lines 21-35)

\(21) **Seed source** This line provides information on the seed sources (premises initially infected to begin simulations). The default of **allFips** will perform a simulation originating in each county present in the premises file (line 11, FLAPS), choosing a random premises within that county as the seed source. Alternatively, you can provide the name of a tab-delimited file containing identifiers (FIPS codes or premisesIDs) from which to seed infection, one line per simulation. You may also seed from multiple premises at once by providing comma-separated premisesIDs in the file (one line per simulation).

\(22) **Seed file type** This describes the information provided in line 21. 

* **fips** = choose 1 premises at random within the FIPS code on a given line in (21) or in all FIPS [default]
* **singlePremises** = use the premisesID provided on a given line in (21)
* **multiplePremises** = use all comma separated premisesIDs on a given line in (21) 

\(23) Unused

\(24) **Susceptibility exponents (mu)** Comma-separated positive numeric values indicating the exponent value for susceptibility for each species/host type, in the same order as listed in line 12. This value is $\mu$ in calculating susceptibility of a premises with *h* animals of a given species as $ah^\mu$ (not including scaling adjustments). Default is 1,1. 

\(25) **Infectiousness exponents (omega)** Comma-separated positive numeric values indicating the exponent value for infectiousness for each species/host type, in the same order as listed in line 12. This value is $\omega$ in calculating infectiousness of a premises with *h* animals of a given species as $bh^\omega$ (not including scaling adjustments). Default is 1,1.

\(26) **Susceptibility constants (a)** Comma-separated positive numeric values indicating the indicating the constant value for susceptibility for each species/host type, in the same order as listed in line 12. This value is *a* in calculating susceptibility of a premises with *h* animals of a given species as $ah^\mu$ (not including scaling adjustments). Default is 1,1.

\(27) **Infectiousness constants (b)** Comma-separated positive numeric values indicating the constant value for infectiousness for each species/host type, in the same order as listed in line 12. This value is *b* in calculating infectiousness of a premises with *h* animals of a given species as $bh^\omega$ (not including scaling adjustments). Default is 5.8,5.8.

\(28) **Kernel type for local (diffusion) spread**

* **0** = Form of kernel describing local spread is $k1/(1 + (distance/k2)^{k3})$ [default]
* **1** = Data-based kernel describing local spread is in file in line 30
* **2** =  Form of kernel describing local spread is $k1/(1+distance/k2)^{k3}$

\(29) **Kernel parameters** Comma-separated positive numeric values for kernel parameters k1, k2, k3 as described in line 28, if applicable. Default is 1.293833e-08,2116.798,2.38

\(30) **Data kernel file** Path and name of file containing data-based kernel referred to in line 28, with columns of distance (in same units as premises x/y) and associated probability of exposure.

\(31) **Latency duration** Two comma-separated positive numeric values describing the mean and variance about a normal distribution for the number of days from premises exposure to infectiousness. Outcomes less than one are set equal to one. Default is 5,0. 

\(32) **Infectiousness duration** Two comma-separated positive numeric values describing the mean and variance about a normal distribution for the number of days from premises infectiousness to immunity. Outcomes less than one are set equal to one. Default is 7,0. 

\(33) Additional functionality under development. Do not change default value [0].

\(34) Additional functionality under development. Do not change default value [*].

\(35) Unused

#### Grid-related settings (lines 36-40)

\(36) **Grid cell file** Filename containing grid cells (will override other options)

\(37) **Grid cell length** Length of cell side for uniform cells (will override density options)
<!-- Does density options mean the ones in line 38? I.e. can i phrase this whole section as the line with the lowest number will over-ride options in cells with higher line numbers? -->

\(38) **Grid cell options** Two comma-separated positive integers specifying a maximum number of premises per cell and minimum side length of a grid cell (same units as premises x/y). This is used to generate grid cells. Default is 500,100000.

\(39) Unused

\(40) Unused

#### Shipment-related settings (lines 41-50)

\(41) **Shipment generation method** Method(s) to generate county-county shipments using USAMM, comma-separated 

* **1** = USAMM version 1 [default] 
* **-1** = shipments off 

\(42) **Shipment method time** Comma-separated timepoints at which to begin using each of the methods above. First timepoint must be 1.

\(43) Additional functionality under development. Do not change default value [0].

\(44) **Parameter files** Output files from USAMM with parameters for a, b, N, s and covariates. Default is inputfiles/usammv1_beef.res, inputfiles/usammv1_dairy.res

\(45) Additional functionality under development. Do not change default value [Q1].

\(46) Additional functionality under development. Do not change default value [1].

\(47) **Origin covariates** Covariates for all origin counties. One file for each species, comma separated. Must have header: FIPS, name1, name2... Names must match name component of covariate parameters in USAMM parameter file (line 44). Default is \*,*.

\(48) **Destination covariates** Covariates for all destination counties. One file for each species, comma separated. Must have header: FIPS, name1, name2... Names must match name component of covariate parameters in USAMM parameter file (line 44). Default is \*,*.

\(49) Disabled

\(50) **Exposed shipments** If enabled, shipping from a farm with "exposed" status will cause the receiver to become exposed as well. 

* **1** = on [default] 
* **0** = off

#### Control-related settings (lines 51-70)

There are two sections defining control actions. Lines 51-61 define attributes of particular control actions. Within this section, attributes for each control action must be listed in the same order as they initially appear in line 51. Lines 63-67 define how and when those control actions are applied.

\(51) **Control types** Comma-separated names for unique control types. Options (which can be combined) are listed below. Based on this line's input there are restricted entries on other lines. These are also listed below. 
<!-- explain the options on a line each? -->

- **\*** turns off control. This option over-rides all others in this section. 

- **shipBan**: Shipment (movement) ban
    - Line 52's section for this control type should be **noLimit**
    - Line 53's section for this control type should be **0**
    - Line 54's section for this control type should be either **county** or **state**
    - Lines 55 and 56's sections for this control type should be **NA**

- **cull**: Cull
    - Line 53's section for this control type will be interpreted as parameters for a daily per premises limit on numbers of animals that can be controlled
    - Line 54's section for this control type should be **premises**

- **vax**: Vaccination
    - Line 53's section for this control type will be interpreted as parameters for a daily per-premises limit on numbers of animals that can be controlled
    - Line 54's section for this control type should be **premises** 

\(52) **Constraint types** Comma-separated names for constraint types to apply to each control type, in the same order as listed in line 51. Options are listed below. Based on this line's input there are restricted entries on other lines. These are also listed below.  

- **noLimit**: Any premises targeted for a particular control action will have that control type immediately implemented in the same timestep
    - Line 53's section for this control type should be **0**
    - Lines 55 and 56's sections for this control type should be **NA**

- **stateSum**: Premises within the same state all draw on a summed pool of resources for that state. Additionally, daily constraints will limit the number of animals per day that can be controlled.
    - Line 53's section for this control type should be the mean and variance for daily per premises limits on numbers of animals that can be controlled 
    - Line 55's section for this control type should be resources to be summed by state
    - Line 56's section for this control type should be **resourceLocs**

- **dailyLimit**: Premises are limited to having some number of animals controlled per day, but there are no additional constraints.
    - Line 53's section for this control type should be the mean and variance for daily per premises limits on numbers of animals that can be controlled

- **maxDistance**: A string combining "maxDistance" with a number denoting the radius for the farthest possible control resource, for example 'maxDistance50000'. This radius is in the same units as premises coordinates.  For a given premises, the nearest available resource within the maximum distance will be assigned. Additionally, daily constraints will limit the number of animals per day that can be controlled. If no resources within the specified distance are available, control will not be performed for that premises
    - Line 53's section for this control type should be the mean and variance for daily per premises limits on numbers of animals that can be controlled
    - Line 55's section for this control type should be resources and their locations. 
    - Line 56's section for this control type should be **resourceLocs**
 
\(53) **Constraint parameters** Semicolon-separated list of numeric parameters to apply to control constraints, in the same order as listed in line 51. Multiple parameters to be applied to a given constraint should be comma-separated. Currently, the only constraint parameters provided here are the mean and variance of a daily per-premises limit on numbers of animals that can be controlled. The constraint parameters provided here are the mean and variance of animals that can be controlled. For example, the values for the stateSum constraint are the mean maximum number of animals per day per premises that can be culled, which is calcualted as 20 cows/hr\*12 hrs/day = 240. The variance is set at 0. The default mean vaccination constraint is 6804, with variance 0.
 
\(54) **Spatial scale** Scale at which each control type should be applied. Possible options are:

- **state**: Currently can only be specified with control type "shipBan"
- **county**: Currently can only be specified with control type "shipBan"
- **premises**: The only valid specification for control types "cull" and "vax"

\(55) **Additional constraint files** Semicolon-separated list of paths and names of files containing additional information to use for control constraints, in the same order as listed in line 51. Multiple files to be applied to a given constraint should be comma-separated. Currently only accepts landfill locations for control type "cull". For control types other than "cull", enter **NA**.

\(56) **Additional constraint type** Types of information provided in line 55, also semicolon separated by type and comma separated for multiple files to be applied to the same control type. Currently only accepts **resourceLocs** for landfill locations for control type "cull". For control types other than "cull", enter **NA**.

\(57) **Mean effectiveness lag** Mean number of days from the control type being implemented to becoming effective (when transmission is affected). Comma-separated positive numeric values, one for each control type in the order listed in line 51, indicating the mean (about a normal distribution). The realized number of days will be drawn from a distribution based on this and line 58, and rounded to the nearest integer. Outcomes less than one are set equal to one. The default for **shipBan** and **cull** is 0, and for **vax** is 11. 

\(58) **Variance effectiveness lag** Variance in the number of days from the control type being implemented to becoming effective (when transmission is affected). Comma-separated positive numeric values, one for each control type in the order listed in line 51, indicating the variance (about a normal distribution). The realized number of days will be drawn from a distribution based on this and line 57, and rounded to the nearest integer. Outcomes less than one are set equal to one. Defaults are 0 for all control types. 

\(59) **Mean inactivation lag** Mean number of days for which the control type is effective before becoming inactive (when transmission is no longer affected). Comma-separated positive numeric values, one for each control type in the order listed in line 51, indicating the mean (about a normal distribution). The realized number of days will be drawn from a distribution based on this and line 60, and rounded to the nearest integer. Outcomes less than one are set equal to one. For control types that should be effectively permanent, set this value higher than the value in line 13. The default for **shipBan** and **cull** is 366, and for **vax** is 183. 

\(60) **Variance inactivation lag** Variance (about a normal distribution) in the number of days for which the control type is effective before becoming inactive (when transmission is no longer affected). Comma-separated positive numeric values, one for each control type in the order listed in line 51, indicating variance (about a normal distribution). The realized number of days will be drawn from a distribution based on this and line 59, and rounded to the nearest integer. Outcomes less than one are set equal to one. Defaults are 0 for all control types.

\(61) **Control effectiveness** Effectiveness (including compliance) of control types. Pairs of proportions (numeric values in [0,1]) separated by commas, semicolon-separated by control type. The first value in the pair is the control type's probability of preventing exposure in a given transmission event. The second value in the pair is the control type's probability of preventing transmission from an infectious premises in a given transmission event. For control type "shipBan", the first and second numbers should be the same and indicate the probability of a a shipment not occurring due to the ban.
<!-- Make sure I understand this, maybe explain better -->

\(62) Unused

\(63) **Trigger type** Triggers that cause a control action to be applied. Options include

- **newPremReportsOverX** = Control will be triggered when the number of premises reported in a timestep exceeds the value specified in line 64. 
- **newRegionReportsOverX** = Control will be triggered when the number of regions (counties or states as specified in line 54) reported in a timestep exceeds the value specified in line 64. 

\(64) **Trigger value** Comma-separated positive integer values, one for each trigger specified in line 63. Trigger thresholds over which control actions in line 63 should be applied.

\(65) **Trigger control** Comma-separated control types to be applied for each trigger specified in line 63 (in the same order). Only control types specified in line 51 can be used. Additional triggers may be added to apply the same control type in different situations (i.e. to cull both reported premises and dangerous contacts).
<!-- Don't get the last sentence -->

\(66) **Response targets** Comma-separated control response targets, in the same order as the triggers specified in line 63. 

* **-1** = apply triggered control type to dangerous contacts of reported farms
* **0** = apply triggered control type to reported farms
* **Positive integer** = radius in units of x/y coordinates of premises (usually meters)

\(67) **Prioritization** Comma-separated prioritization methods for adding to control waitlists, in the same order as the triggers specified in line 63. Currently only accepts **earliest** for each control type.

\(68) Unused

\(69) Unused 

\(70) Unused


#### Reporting and Tracing (DC) settings (lines 71-75)

\(71) **Index reporting lag** Two comma-separated positive numeric values describing the mean and variance about a normal distribution for the number of days from index premises exposure (the first case) to reporting. Outcomes less than one are set equal to one. Default is 15,0. 

\(72) **Reporting lag** Two comma-separated positive numeric values describing the mean and variance about a normal distribution for the number of days from premises exposure to reporting, if the premises has NOT been designated as a dangerous contact. Outcomes less than one are set equal to one. Default is 8,0.

\(73) **DC reporting lag** Two comma-separated positive numeric values describing the mean and variance about a normal distribution for the number of days from premises exposure to reporting, if the premises HAS been designated as a dangerous contact. Outcomes less than one are set equal to one. Default is 2,0.

\(74) **DC scaling** Dangerous contact scaling parameters relative to risk. Only used if at least one value in line 66 is -1, i.e. if dangerous contacts are a specified control target. Comma-separated pairs of a string and a numeric value >=1, semicolon separated. The string is the disease status for which the scaling constant will apply, and the number is the scaling constant. Generally, there will be two pairs, one each for statuses **sus** and **exp**. The value for 'exp' generally will be higher than the value for 'sus', indicating that premises are more likely to be identified as dangerous contacts if they truly are exposed at the time of this identification than if they are not. Default is **sus,4; exp,5**

\(75) Unused 
