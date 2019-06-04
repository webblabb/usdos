#include "File_manager.h"

File_manager::File_manager()
{
	verbose = verboseLevel;
}

File_manager::~File_manager()
{
	delete params.kernel;
}

/// Returns a string to output to run log file, containing batchDateTime (provided as argument), and contents of pv (config file contents)
const std::string File_manager::getSettings(std::string& bdt)
{
	std::string output = bdt;
	for (int i=1; i<=75; i++){
		output += "\t";
		output += pv[i];
	}
	output += "\n";
	return output;
}

/// Reads configuration file, stores parameter values in a character vector.
/// Checks validity of values and groups closely related parameters.
void File_manager::readConfig(std::string& cfile)
{
	pv.emplace_back("0"); // fills in [0] so that line numbers match up with elements
	// Read in file and store in parameter vector pv (private class variable)
	std::ifstream f(cfile);
	if(f.is_open()){
		while(!f.eof()){ // until end of file
			std::string line;
			std::stringstream line2;
			// First get a whole line from the config file
			std::getline(f, line);
			// Put that in a stringstream (required for getline) and use getline again
			// to only read up until the comment character (delimiter).
			line2 << line;
			std::getline(line2, line, '#');
			// If there is whitespace between the value of interest and the comment character
			// this will be read as part of the value. Therefore strip it of whitespace first.
			line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
			if(line.size() != 0){pv.emplace_back(line);}
		}
		if (pv.size()!=76){std::cout<<"Warning: expected configuration file with 75 lines, loaded file with "<<pv.size()-1<<" lines."<<std::endl;
		} else {
		if (verbose >0){
		std::cout << "Configuration file loaded with "<<pv.size()-1<<" lines."<<std::endl;}
		}
		// Check for consistencies/requirements in parameter vectors, spit out warnings/errors and exit if needed
		bool exitflag = 0;
		bool checkExit = 0;

		// Batch name
		params.batch = pv[1];
		// Outputs on/off
		if (pv[2]=="*" || pv[3]=="*"|| pv[4]=="*"|| pv[5]=="*"|| pv[6]=="*"){
			std::cout << "Warning: (config 1-5): Output should be specified - setting unspecified (*) to off." << std::endl;
			if (pv[2]=="*"){pv[2]="0";}
			if (pv[3]=="*"){pv[3]="0";}
			if (pv[4]=="*"){pv[4]="0";}
			if (pv[5]=="*"){pv[5]="0";}
			if (pv[6]=="*"){pv[6]="0";}
		}
		params.printSummary = stringToNum<int>(pv[2]);
		params.printDetail = stringToNum<int>(pv[3]);
		params.printCells = stringToNum<int>(pv[4]);
		params.printShipments = stringToNum<int>(pv[5]);
		params.printControl = stringToNum<int>(pv[6]);
		// pv[7] ... pv[10]
		// Premises file
		if (pv[11]=="*"){
			std::cout << "ERROR (config 11): No premises file specified." << std::endl; exitflag=1;}
		params.premFile = pv[11];
		// Species
		if (pv[12]=="*"){
			std::cout << "ERROR (config 12): No species list provided." << std::endl; exitflag=1;}
		params.species = stringToStringVec(pv[12]);
		// Timesteps
		params.timesteps = stringToNum<int>(pv[13]);
		if (params.timesteps<1){std::cout << "Warning (config 13): Number of timesteps must be 1 or more. Setting number of timesteps to 365." << std::endl;
			params.timesteps = 365;}
		// Max number of exposed premises
		params.useMaxPrems = 0;
		if (pv[14] != "*"){
			params.maxInfectiousPrems = stringToNum<int>(pv[14]);
			if (params.maxInfectiousPrems<1){
				std::cout << "Warning (config 14): Number of maximum premises must be 1 or more. Value will not be set." << std::endl;
			} else if (params.maxInfectiousPrems>=1){
				params.useMaxPrems = 1;
			}
		}
		// Verbose level
		params.verboseLevel = stringToNum<int>(pv[15]);
		if (params.verboseLevel!=0 && params.verboseLevel!=1 && params.verboseLevel!=2){
			std::cout << "Warning (config 15): Verbose option must be 0, 1 or 2. Setting option to off." << std::endl;
			params.verboseLevel = 0;}
		verbose = params.verboseLevel;
		// Pairwise on
		params.pairwiseOn = stringToNum<int>(pv[16]) ;
		if (params.pairwiseOn!=0 && params.pairwiseOn!=1){
			std::cout << "ERROR (config 16): Pairwise algorithm must be 1 (on) or 0 (off)." << std::endl; exitflag=1;}
		// Reverse x/y
		params.reverseXY = stringToNum<int>(pv[17]);
		if (params.reverseXY!=0 && params.reverseXY!=1){
			std::cout << "ERROR (config 17): Reversing x/y must be 0 (y first) or 1 (x first)." << std::endl; exitflag=1;}
		//County data file
		if (pv[18] == "*"){
            std::cout << "ERROR (config 18): No county data file specified." << std::endl; exitflag=1;}
		params.fipsFile = pv[18];
		//Day of the year to start simulation/generation.
        if(pv[19]=="*"){
            std::cout << "ERROR (config 19): No start day option given in config file." << std::endl;}
        params.start_day_option = stringToNum<int>(pv[19]);
        if(params.start_day_option == 0) {params.start_day = rand_int(1, 365);}
        else {params.start_day = params.start_day_option;}
        //pv[20];
		// Infectious seed file
		if (pv[21]=="*"){
			std::cout << "ERROR (config 21): No infectious premises seed source specified." << std::endl; exitflag=1;}
		params.seedSource = pv[21];
		// Infection seed method
		params.seedSourceType = pv[22];
		if (params.seedSourceType.compare("fips") == 0 &&
				params.seedSourceType.compare("singlePremises") == 0 &&
				params.seedSourceType.compare("multiplePremises") == 0){
		std::cout << "ERROR (config 22): Seed source type must be 'fips','singlePremises', or 'multiplePremises'." << std::endl; exitflag=1;}
		// pv[23] unused
		// Susceptibility exponents by species
		std::vector<double> tempVec1 = stringToNumVec(pv[24]);
		checkExit = checkPositive(tempVec1, 24); if (checkExit==1){exitflag=1;}
		if (tempVec1.size() != params.species.size()){
			std::cout<<"ERROR (config 12 & 24): Different numbers of species and susceptibility exponents provided: "<<params.species.size()<<" species and "
			<<tempVec1.size()<<" exponents." <<std::endl; exitflag=1;}
		// Infectiousness exponents by species
		std::vector<double> tempVec2 = stringToNumVec(pv[25]);
		checkExit = checkPositive(tempVec2, 25); if (checkExit==1){exitflag=1;}
		if (tempVec2.size() != params.species.size()){
			std::cout<<"ERROR (config 25): Different numbers of species and infectiousness exponents provided." <<std::endl; exitflag=1;}
		// Susceptibility constants by species
		std::vector<double> tempVec3 = stringToNumVec(pv[26]);
		checkExit = checkPositive(tempVec3, 26); if (checkExit==1){exitflag=1;}
		if (tempVec3.size() != params.species.size()){
			std::cout<<"ERROR (config 12 & 26): Different numbers of species and susceptibility constants provided: "<<params.species.size()<<" species and "
			<<tempVec3.size()<<" constants." <<std::endl; exitflag=1;}
		// Infectiousness constants by species
		std::vector<double> tempVec4 = stringToNumVec(pv[27]);
		checkExit = checkPositive(tempVec4, 27); if (checkExit==1){exitflag=1;}
		if (tempVec4.size() != params.species.size()){
			std::cout<<"ERROR (config 12 & 27): Different numbers of species and infectiousness constants provided: "<<params.species.size()<<" species and "
			<<tempVec4.size()<<" constants." <<std::endl; exitflag=1;}
		// Map exponents and constant values to species
		int i = 0;
		for (auto &sp:params.species){
			params.susExponents[sp] = tempVec1.at(i);
			params.infExponents[sp] = tempVec2.at(i);
			params.susConsts[sp] = tempVec3.at(i);
			params.infConsts[sp] = tempVec4.at(i);
			++i;
		}
		// Kernel type & parameters
		if (pv[28]!="0" && pv[28]!="1" && pv[28]!="2"){std::cout << "ERROR (config 28): Kernel type must be 0 (power law), 1 (data-based), or 2 (kernel4)." << std::endl; exitflag=1;}
		params.kernelType = stringToNum<int>(pv[28]);
		params.kernelParams = stringToNumVec(pv[29]);
		if ((params.kernelType == 0 || params.kernelType == 2) && params.kernelParams.size() < 3){
			std::cout<<"ERROR (config 29): For kernel type 0 or 2, three kernel parameters are required." <<std::endl; exitflag=1;}
		// Data-based kernel file (checked against kernel option)
		params.dataKernelFile = pv[30];
		if (params.dataKernelFile!="*" && params.kernelType!=1){std::cout<<"Warning: Data kernel file provided but kernel type is not set to 1; file will be ignored." << std::endl;}
		if (params.kernelType == 1 && params.dataKernelFile=="*"){std::cout<<"ERROR (config 28): Kernel type 1 (data-based) requires file in config 30."<< std::endl; exitflag=1;}
		// Latency parameters
		checkExit = checkMeanVar(pv[31],31,"latency"); if (checkExit==1){exitflag=1;} // if exit triggered by this check, set exitflag=1
		std::vector<double> tempVec = stringToNumVec(pv[31]);
		params.latencyParams = std::make_tuple(tempVec[0],tempVec[1]);
		// Infectiousness parameters
		checkExit = checkMeanVar(pv[32],32,"infectiousness"); if (checkExit==1){exitflag=1;} // if exit triggered by this check, set exitflag=1
		tempVec = stringToNumVec(pv[32]);
		params.infectiousParams = std::make_tuple(tempVec[0],tempVec[1]);
		// Partial transition off (0), FMD like (1) and bTB like (2)
        if (pv[33]!="0" && pv[33]!="1" && pv[33]!="2"){std::cout << "ERROR (config 33): Partial transition flag must be 0 (off), 1 (FMD like), or 2 (bTB like)." << std::endl; exitflag=1;}
        params.partial=stringToNum<int>(pv[33]);

		if (pv[33]=="2"){
			std::cout << "ERROR (config 33): bTB partial transition not coded yet." << std::endl; exitflag=1;}

        // Partial transition parameter values
        if (pv[33]=="1"||pv[33]=="2"){
         params.partialParams = stringToNumVec(pv[34]);
            checkExit = checkPositive(params.partialParams, 34); if (checkExit==1){exitflag=1;}
            if (params.partialParams.size() != 6){
                std::cout<<"ERROR (34): Expecting 6 parameter values for partial transition extension." <<std::endl; exitflag=1;}
            if(std::get<1>(params.latencyParams)!=0){
                 std::cout<<"ERROR (31): Variance of days from premises exposure to infectiousness must be 0 if partial transistion flag is on." <<std::endl; exitflag=1;
                }
        }
        if (pv[33]=="2"){std::cout << "ERROR (config 33): Coming soon to a command line near you: bTB like infection for US-DOS." << std::endl; exitflag=1;}

//        //parameters related to open population process
//        if(pv[33]=="2"&pv[35]=="*"){
//            std::cout<<"ERROR (35): Missing aruguments, expecting population paramaters for bTB like infection." <<std::endl; exitflag=1;
//            }
//
//        if(pv[33]=="2"){
//            std::vector<double> tempVec5 = stringToNumVec(pv[35]);
//            checkExit = checkPositive(tempVec5, 35); if (checkExit==1){exitflag=1;}
//            int j = 0;
//            for (auto &sp:params.species){
//                params.popParams[sp] = tempVec5.at(j);
//                ++j;
//            }
//        }

		// Grid settings
		if (pv[36]=="*" && pv[37]=="*" && pv[38]=="*"){std::cout << "ERROR (config 36-38): No grid cell parameters specified." << std::endl; exitflag=1;}
		params.cellFile = pv[36];
		params.uniformSide = stringToNum<double>(pv[37]);
		params.densityParams = stringToIntVec(pv[38]);
			checkExit = checkPositive(params.densityParams, 38); if (checkExit==1){exitflag=1;}
			if ((params.densityParams).size()!=2){std::cout << "ERROR (config 38): Two parameters required for grid creation by density." << std::endl; exitflag=1;}
		//pv[39] ... pv[40]

		// Shipping methods and times
		if (pv[41]=="*"){std::cout << "ERROR (config 41): No county-level shipment method(s) specified." << std::endl; exitflag=1;}
        params.shipMethods = std::stoi(pv[41]);
        params.shipments_on = 0;
        params.shipment_kernel = "off";
        params.usamm_version = 0;
        if(params.shipMethods > 0 ){ //The following options are only of interest if shipments are not turned off.
            params.shipments_on = 1;
            switch(params.shipMethods)
            {
            case 1:
                params.shipment_kernel = "original USAMM";
                params.usamm_version = 1;
                break;
            case 2:
                params.shipment_kernel = "original USAMM";
                params.usamm_version = 2;
                break;
            case 3:
                params.shipment_kernel = "power_exp";
                params.usamm_version = 2;
                break;
            case 4:
                params.shipment_kernel = "one_minus_exp";
                params.usamm_version = 2;
                break;
            case 5:
                params.shipment_kernel = "local";
                params.usamm_version = 2;
                break;
            default:
                params.shipment_kernel = "::UNKNOWN, CHECK CONFIG FILE::";
                params.usamm_version = -1;
                break;
            }
            if (pv[42]=="*"){std::cout << "ERROR (config 42): No shipment method start time(s) specified." << std::endl; exitflag=1;}
            params.shipMethodTimeStarts = stringToIntVec(pv[42]);
//            if (params.shipMethods.size() != params.shipMethodTimeStarts.size()){
//                std::cout << "ERROR (config 41-42): Number of methods and start times provided must be the same ("<<params.shipMethods.size()<<" method(s) and "<<params.shipMethodTimeStarts.size()<<" start time(s) provided)."<<std::endl; exitflag=1;}
            if (params.shipMethodTimeStarts.at(0) != 1){std::cout << "ERROR (config 42): First county-level shipment method start time must be 1."<<std::endl; exitflag=1;}
             // check that each value is between 1 and number of timesteps
            bool rangeflag = 0;
            for (auto& t:params.shipMethodTimeStarts){if(t<1 || t>params.timesteps){rangeflag=1;}}
            if (rangeflag){std::cout << "Warning (config 42): Simulation timespan (config 13) is shorter than largest timepoint for county-level shipment methods - some methods may not be used." << std::endl;}
            params.shipMethodTimeStarts.emplace_back(params.timesteps+1); // add this b/c whichElement function needs a maximum (element won't be used)
            if (pv[43]=="*"){std::cout << "ERROR (config 43): No premises-level shipment assignment method specified." << std::endl; exitflag=1;}
            params.shipPremAssignment = stringToNum<int>(pv[43]);
//            if(params.shipMethods.front() == 2 and pv[44]=="*"){std::cout << "ERROR (config 44): No USAMM parameter file(s) specified." << std::endl; exitflag=1;}
            params.USAMM_parameter_files = stringToStringVec(pv[44]);
            if(pv[44]!="*" and params.USAMM_parameter_files.size() != params.species.size()){std::cout << "ERROR (config 12 & 44): The number of USAMM parameter files and number of species in the config file does not match." << std::endl; exitflag=1;}
//            if(params.shipMethods.front() == 2 and pv[45]=="*"){std::cout << "ERROR (config 45): No USAMM temporal ordering specified." << std::endl; exitflag=1;}
            params.USAMM_temporal_order = stringToStringVec(pv[45]);
//            if(params.shipMethods.front() == 2 and pv[46]=="*"){std::cout << "ERROR (config 46): No USAMM temporal start points specified." << std::endl; exitflag=1;}
            params.USAMM_temporal_start_points = stringToIntVec(pv[46]);
            params.USAMM_temporal_n_days.resize(params.USAMM_temporal_start_points.size(), 0);
            int counter = 366;
            for(int i = params.USAMM_temporal_start_points.size()-1; i >= 0; i--)
            {
                params.USAMM_temporal_n_days[i] = counter-params.USAMM_temporal_start_points[i];
                counter = params.USAMM_temporal_start_points[i];
            }

//            if(params.shipMethods.front() == 2 and pv[47]=="*"){std::cout << "ERROR (config 47): No USAMM origin covariate file(s) specified." << std::endl; exitflag=1;}
            params.USAMM_ocov_files = stringToStringVec(pv[47]);
//            if(params.shipMethods.front() == 2 and pv[48]=="*"){std::cout << "ERROR (config 48): No USAMM destination covariate file(s) specified." << std::endl; exitflag=1;}
            params.USAMM_dcov_files = stringToStringVec(pv[48]);
//            if(params.shipMethods.front() == 2 and pv[49]=="*"){std::cout << "ERROR (config 49): No USAMM supernode files specified." << std::endl; exitflag=1;}
            params.USAMM_supernode_files = stringToStringVec(pv[49]);

            if (params.USAMM_ocov_files.size() != params.species.size()){
                std::cout<<"ERROR (config 12 & 47): Different numbers of origin covariate input files: "<<params.species.size()<<" species and "
                         <<params.USAMM_ocov_files.size()<<" origin covariate files." <<std::endl; exitflag=1;}
            if (params.USAMM_dcov_files.size() != params.species.size()){
                std::cout<<"ERROR (config 12 & 47): Different numbers of destination covariate input files: "<<params.species.size()<<" species and "
                         <<params.USAMM_dcov_files.size()<<" destination covariate files." <<std::endl; exitflag=1;}
            if (params.USAMM_supernode_files.size() != params.species.size()){
                std::cout<<"ERROR (config 12 & 47): Different numbers of supernode input files: "<<params.species.size()<<" species and "
                         <<params.USAMM_supernode_files.size()<<" supernode files." <<std::endl; exitflag=1;}


			params.exposed_shipments = false;
			params.statuses_to_generate_shipments_from = {"inf"};
			std::string exposed_shipment_option = pv[50];
			if(exposed_shipment_option.compare("1") == 0)
            {
                params.exposed_shipments = true;
                params.statuses_to_generate_shipments_from.push_back("exp");
            }
        }

		// Control - type names
		params.control_on = true;
		if (pv[51]=="*" || pv[63]=="*"){
		params.control_on = false;}
		params.dangerousContacts_on = 0;

if(params.control_on == true){
		params.controlTypes = stringToStringVec(pv[51]);		
			std::string options_ct = "cull,vax,shipBan";
			if (verbose>1){std::cout<< "control type options: " << pv[51]<<std::endl;}
			
			std::vector<std::string> options = stringToStringVec(options_ct);			
			checkExit = limitedValues(params.controlTypes, options, 51); if (checkExit==1){exitflag=1;}
			
		// Control - constraint function names
		std::vector<std::string> cFuncs = stringToStringVec(pv[52]);
		if (verbose>1){std::cout<< "constraint function names options: " << pv[52]<<std::endl;}
		
			if (cFuncs.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 51 and 52): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
			std::string options_cf = "noLimit,stateSum,dailyLimit,nationalLimit";
			options = stringToStringVec(options_cf);
			checkExit = limitedValues(cFuncs, options, 52); if (checkExit==1){exitflag=1;}
		// Control - constraint function parameters
		std::vector<std::string> cFuncParams = semicolonStringToStringVec(pv[53]);
			if (cFuncParams.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 51 and 53): Different number of arguments and control types. Provide 0 for constraintFunctions = 'noLimit'. Or are parameter sets not separated by semicolon?" << std::endl;
			exitflag=1;
			}
			if (verbose>1){std::cout<< "constraint paramters options: " << pv[53]<<std::endl;}
			
		// Control - spatial scale for each type
		std::vector<std::string> spatialScales = stringToStringVec(pv[54]);
			std::string options_scales = "premises,county,state";
			options = stringToStringVec(options_scales);
			checkExit = limitedValues(spatialScales, options, 54); if (checkExit==1){exitflag=1;}
			if (spatialScales.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 54): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
			if (verbose>1){std::cout<< "spatial scale options: " << pv[54]<<std::endl;}
			
		// Control - external files relating to constraint functions
		std::vector<std::string> cFuncFiles = semicolonStringToStringVec(pv[55]);
		if (verbose>1){std::cout<< "control files for constraint options: " << pv[55]<<std::endl;}
		
			if (cFuncFiles.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 55): Different number of arguments and control types. Are sets not separated by semicolon?" << std::endl;
			exitflag=1;
			}
		// Control - types for each external file relating to constraint functions
		std::vector<std::string> cFuncFileTypes = semicolonStringToStringVec(pv[56]);
		if (verbose>1){std::cout<< "control files type options: " << pv[56]<<std::endl;}
			std::string options_filetypes = "NA,resourceLocs,resourceBoosts"; //  added time or outbreak dependent supply boosts
			options = stringToStringVec(options_filetypes);
			checkExit = limitedValues(cFuncFileTypes, options, 56); if (checkExit==1){exitflag=1;}
			if (cFuncFileTypes.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 56): Different number of arguments and control types. Are sets not separated by semicolon?" << std::endl;
			exitflag=1;
			}
		// Control - time lags for each control: initiation to effective
		std::vector<double> implementToEffectiveMeans = stringToNumVec(pv[57]);
			checkExit = checkPositive(implementToEffectiveMeans, 57); if (checkExit==1){exitflag=1;}
			if (implementToEffectiveMeans.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 57): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
		// Control - time lags for each control: initiation to effective
		std::vector<double> implementToEffectiveVars = stringToNumVec(pv[58]);
			checkExit = checkPositive(implementToEffectiveVars, 58); if (checkExit==1){exitflag=1;}
			if (implementToEffectiveVars.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 58): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
		// Control - time lags for each control: effective to inactive
		std::vector<double> effectiveToInactiveMeans = stringToNumVec(pv[59]);
			checkExit = checkPositive(effectiveToInactiveMeans, 59); if (checkExit==1){exitflag=1;}
			if (effectiveToInactiveMeans.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 59): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
		// Control - time lags for each control: effective to inactive
		std::vector<double> effectiveToInactiveVars = stringToNumVec(pv[60]);
			checkExit = checkPositive(effectiveToInactiveVars, 60); if (checkExit==1){exitflag=1;}
			if (effectiveToInactiveVars.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 60): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
		// Control - effectiveness of each control, including compliance
		std::vector<std::string> effAll = semicolonStringToStringVec(pv[61]);
		if (effAll.size() != params.controlTypes.size()){
			std::cout << "ERROR (config 53 and 61): Different number of arguments and control types." << std::endl;
			exitflag=1;
			}
		// pv[62]
		// Control rules
		std::vector<std::string> triggers = stringToStringVec(pv[63]);
			std::string options_triggers = "newPremReportsOverX,newRegionReportsOverX";
			options = stringToStringVec(options_triggers);
			checkExit = limitedValues(triggers, options, 63); if (checkExit==1){exitflag=1;}
		std::vector<double> thresholds = stringToNumVec(pv[64]);
			checkExit = checkPositive(thresholds, 64); if (checkExit==1){exitflag=1;}
			if (thresholds.size() != triggers.size()){
			std::cout << "ERROR (config 63 and 64): Different number of arguments and control triggers." << std::endl;
			exitflag=1;
			}
		std::vector<std::string> actions = stringToStringVec(pv[65]);
			options = stringToStringVec(options_ct); // control types from line 51
			checkExit = limitedValues(actions, options, 65); if (checkExit==1){exitflag=1;}
			if (actions.size() != triggers.size()){
			std::cout << "ERROR (config 63 and 65): Different number of arguments and control triggers." << std::endl;
			exitflag=1;
			}
		std::vector<double> targets = stringToNumVec(pv[66]);	// should be -1, 0, or positive
			if (targets.size() != triggers.size()){
			std::cout << "ERROR (config 63 and 66): Different number of arguments and control triggers." << std::endl;
			exitflag=1;
			}
			for (auto&t:targets){
				if (t==-1){params.dangerousContacts_on = 1;}
			}

// Should also check that any regional control types (scale != premises) can only have target=0
		std::vector<std::string> orders = stringToStringVec(pv[67]);
			std::string options_orders = "earliest"; // can add "largest, closest" later
			options = stringToStringVec(options_orders);
			checkExit = limitedValues(orders, options, 51); if (checkExit==1){exitflag=1;}
			if (orders.size() != triggers.size()){
			std::cout << "ERROR (config 63 and 67): Different number of arguments and control triggers." << std::endl;
			exitflag=1;
			}

		//pv[68]... pv[70]

		// Reporting lags - index
		checkExit = checkMeanVar(pv[71],71,"index reporting"); if (checkExit==1){exitflag=1;}
		tempVec = stringToNumVec(pv[71]);
		params.indexReportLag = std::make_tuple(tempVec[0],tempVec[1]);
		// Reporting lags - normal, not dangerous contact
		checkExit = checkMeanVar(pv[72],72,"non-DC reporting"); if (checkExit==1){exitflag=1;}
		tempVec = stringToNumVec(pv[72]);
		params.nonDCReportLag = std::make_tuple(tempVec[0],tempVec[1]);
		// Reporting lags - dangerous contacts
		checkExit = checkMeanVar(pv[73],73,"DC reporting"); if (checkExit==1){exitflag=1;}
		tempVec = stringToNumVec(pv[73]);
		params.dcReportLag = std::make_tuple(tempVec[0],tempVec[1]);
		// Dangerous contact scaling
		std::vector<std::string> dcScaling = semicolonStringToStringVec(pv[74]);
		if (params.dangerousContacts_on == 1){
			params.maxDCScale = 0;
			for (auto& d:dcScaling){ // d is std::string that looks like "sus,1.3"
				std::vector<std::string> temp = stringToStringVec(d);
				if (temp.size() != 2){
					std::cout<<"ERROR (config 74): Only two arguments (status and scale) per semicolon-separated set. Exiting...";
					exitflag=1;
				} else { // has 2 elements, check that first is valid status
				  std::string statusAtSourceReport = temp.at(0);
					std::vector<std::string> statusAtSourceReportVec; statusAtSourceReportVec.emplace_back(statusAtSourceReport);
					std::string options_dc = "sus,exp,inf";
					options = stringToStringVec(options_dc);
					checkExit = limitedValues(statusAtSourceReportVec, options, 74); if (checkExit==1){exitflag=1;}
					// check that second is positive number
					double dcScale = 0;
					try{
						dcScale = stringToNum<double>(temp.at(1));
					} catch(std::exception& e){
						std::cout<<"ERROR (config 74): Second argument in set must be numeric. Exiting...";
						exitflag=1;
					}
					if (dcScale<0){
						std::cout<<"ERROR (config 74): Second argument in set must be >0. Exiting...";
						exitflag=1;
					}
					params.dcRiskScale[statusAtSourceReport] = dcScale;
					// save largest of scale values
					if (dcScale > params.maxDCScale){
						params.maxDCScale = dcScale;
					}
 				} // end "if 2 arguments provided in set"
			} // end "for each set of arguments"
		} // end "if dangerousContacts_on"

		// Put control parameters into proper containers
		for (int ct2 = 0; ct2 < params.controlTypes.size(); ++ct2){
			std::string c_type = params.controlTypes.at(ct2);

			params.controlScales[c_type] = spatialScales.at(ct2);
			params.constraintFunctions[c_type] = cFuncs.at(ct2);
			params.constraintFuncParams[c_type] = stringToNumVec(cFuncParams.at(ct2));
			params.constraintFuncFiles[c_type] = stringToStringVec(cFuncFiles.at(ct2));
			params.constraintFuncFileTypes[c_type] = stringToStringVec(cFuncFileTypes.at(ct2));
			params.implementToEffectiveLag[c_type] = std::make_tuple(implementToEffectiveMeans.at(ct2),
				implementToEffectiveVars.at(ct2));
			params.effectiveToInactiveLag[c_type] = std::make_tuple(effectiveToInactiveMeans.at(ct2),
				effectiveToInactiveVars.at(ct2));
			std::vector<double> effectVector = stringToNumVec(effAll.at(ct2));
			checkExit = checkZeroToOne(effectVector, 61); if (checkExit==1){exitflag=1;}
			params.effectiveness[c_type] = std::make_tuple(effectVector.at(0), effectVector.at(1));
		}	// end "for each control type"

		// Put control rules into containers
		for (int cr = 0; cr < triggers.size(); cr++){
			params.controlRules.emplace_back(
				controlRule{
					triggers[cr],
					thresholds[cr],
					actions[cr],
					targets[cr],
					orders[cr],
					targets[cr]*targets[cr] // radius squared, stored pre-calculated for speed
				});
			if (targets[cr] == -1){
				params.dcControlTypes.emplace_back(actions[cr]);
			}
		}
} // end "if control_on==true"

if (exitflag){
	std::cout << "Exiting..." << std::endl;
	exit(EXIT_FAILURE);
}
		// Construct kernel object
		switch (params.kernelType)
		{
			case 0:{
				params.kernel = new Local_spread(0, params.kernelParams);
				break;
			}
			case 1:{
				params.kernel = new Local_spread(1, params.dataKernelFile);
				break;
			}
			case 2:{
				params.kernel = new Local_spread(2, params.kernelParams);
				break;
			}
		}

if (verbose>0){std::cout<<"Parameter loading complete."<<std::endl;}


	} else { // if file not found
		std::cout << "ERROR: Configuration file not found: " << cfile << std::endl <<
		"Exiting..." << std::endl;
		exit(EXIT_FAILURE);
	}
}

/// Use with any argument meant to be in the form (mean,variance).
/// Checks that string s contains two positive arguments.
/// If only one argument is provided, that will be used as the mean (with warning).
/// Outputs error information using line number (lineNum) and parameter description (paramDesc).
/// Returns exitFlag true/false.
bool File_manager::checkMeanVar(std::string& s, int lineNum, std::string paramDesc)
{
	bool exitflag = 0;
	std::vector<double> tempVec;
	tempVec = stringToNumVec(s);
		if (tempVec.size()==0){std::cout << "ERROR (config "<<lineNum<<"): No "<<paramDesc<<" parameters specified." << std::endl; exitflag=1;
		} else if (tempVec.size()==1){std::cout << "Warning (config "<<lineNum<<"): Only one "<<paramDesc<<" parameter provided, mean will be set to "<<tempVec[0]<<" and variance to 0."<<std::endl;
			s.append(",0");
		} else if (tempVec.size()>2){std::cout << "ERROR (config "<<lineNum<<"): Only 2 "<<paramDesc<<" parameters (mean,variance) permitted, "<<tempVec.size()<<" were provided." << std::endl; exitflag=1;}

		bool checkPos = checkPositive(tempVec, lineNum); if (checkPos==1){exitflag=1;}
	return exitflag;
}

bool File_manager::checkPositive(std::vector<int>& tempVec, int lineNum)
{
	bool exitflag = 0;
	unsigned int it = 0;
	while (it < tempVec.size() && exitflag ==0){
		auto tv = tempVec[it];
		if(tv<0){
			std::cout << "ERROR (config "<<lineNum<<"): All parameters must be positive. " << std::endl;
			exitflag=1;
		}
		it++;
	}
	return exitflag;
}

// overloaded for doubles
bool File_manager::checkPositive(std::vector<double>& tempVec, int lineNum)
{
	bool exitflag = 0;
	unsigned int it = 0;
	while (it < tempVec.size() && exitflag ==0){
		auto tv = tempVec[it];
		if(tv<0){
			std::cout << "ERROR (config "<<lineNum<<"): All parameters must be positive. " << std::endl;
			exitflag=1;
		}
		it++;
	}
	return exitflag;
}

/// Check that value is between 0 and 1
bool File_manager::checkZeroToOne(std::vector<double>& tempVec, int lineNum)
{
	bool exitflag = 0;
	unsigned int it = 0;
	while (it < tempVec.size() && exitflag ==0){
		auto tv = tempVec[it];
		if(tv<0 || tv>1){
			std::cout << "ERROR (config "<<lineNum<<"): All values must be between 0 and 1. " << std::endl;
			exitflag=1;
		}
		it++;
	}
	return exitflag;
}

/// Use with entries meant to be a list of limited values:
/// Checks that string s contains two positive arguments
/// Outputs error information using line number (lineNum)
/// Returns exitFlag true/false
bool File_manager::limitedValues(std::vector<std::string>& svec,
	std::vector<std::string>& options, int lineNum)
{
	bool exitflag = 0;
	std::vector<std::string>::iterator it = svec.begin();

	while (exitflag == 0 && it != svec.end()){
		bool hasMatch = 0;
		std::vector<std::string>::iterator it2 = options.begin();
		while (hasMatch == 0 && it2 != options.end()){
			if ((*it) == (*it2)){hasMatch = 1;}
			it2++;
		}
		if (hasMatch == 0){exitflag = 1;}
		it++;
	}

	if (exitflag == 1){
		std::cout << "ERROR (config "<<lineNum<<"): All values must be one of the following";
		for (auto& a:options){
			std::cout << ", " << a;
		}
		std::cout << std::endl;
	}

	return exitflag;
}
