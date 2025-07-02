//© 2024 Colorado State University
#include "File_manager.h"
#include <fstream>

const std::map<DiseaseStatus, std::string> Parameters::DiseaseStatusToString =
    { {DiseaseStatus::SUS, "sus"}, {DiseaseStatus::EXP, "exp"},
      {DiseaseStatus::INF, "inf"}, {DiseaseStatus::IMM, "imm"},
      {DiseaseStatus::BTB, "btb"}, {DiseaseStatus::REC, "rec"} };
const std::map<std::string, DiseaseStatus> Parameters::StringToDiseaseStatus =
    { {"sus", DiseaseStatus::SUS}, {"exp", DiseaseStatus::EXP},
      {"inf", DiseaseStatus::INF}, {"imm", DiseaseStatus::IMM},
      {"btb", DiseaseStatus::BTB}, {"rec", DiseaseStatus::REC} };

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
	for (size_t i=1; i<pv.size(); i++){
		output += "\t";
		output += pv[i];
	}
	output += "\n";
	return output;
}

/// Reads configuration file, stores parameter values in a character vector.
/// Checks validity of values and groups closely related parameters.
void File_manager::readConfig(std::string& cfile, bool generate_network)
{
    params.generate_shipment_network = generate_network;
	pv.emplace_back("0"); // fills in [0] so that line numbers match up with elements
	// Read in file and store in parameter vector pv (private class variable)
	std::ifstream f(cfile);
	if(f.is_open()){
		while(!f.eof()){ // until end of file
			std::string line = "";
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
		if (pv.size()!=150){std::cout<<"Warning: expected configuration file with 149 lines, loaded file with "<<pv.size()-1<<" lines."<<std::endl;
		} else {
		if (verbose >0){
		std::cout << "Configuration file loaded with "<<pv.size()-1<<" lines."<<std::endl;}
		}
		// Check for consistencies/requirements in parameter vectors, spit out warnings/errors and exit if needed
		bool exitflag = 0;
		bool checkExit = 0;

		// Batch name
		params.batch = pv[1];
        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer,80,"_%Y%b%d_%H%M",timeinfo);
        std::string str(buffer);
        params.batchDateTime = params.batch + str;

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
		params.printDiagnostic = stringToNum<int>(pv[7]);
		// pv[8] ... pv[10]
		// Premises file
		if (pv[11]=="*"){
			std::cout << "ERROR (config 11): No premises file specified." << std::endl; exitflag=1;}
		params.premFile = pv[11];
		// Species
		if (pv[12]=="*"){
			std::cout << "ERROR (config 12): No species list provided." << std::endl; exitflag=1;}
		params.species = stringToStringVec(pv[12]);
		for(const auto& sp_name : params.species)
        {
            if(sp_name != "beef" and sp_name != "dairy" and sp_name != "swine")
            {
                std::cout << "ERROR (config 12): List of species must consist of \"beef\", \"dairy\" or \"swine\". Found " << sp_name << "." << std::endl; exitflag=1;
            }
        }

		// Timesteps
		if(generate_network)
        {
            params.timesteps = 365;
        }
        else
        {
            params.timesteps = stringToNum<int>(pv[13]);
            if (params.timesteps<1){std::cout << "Warning (config 13): Number of timesteps must be 1 or more. Setting number of timesteps to 365." << std::endl;
                params.timesteps = 365;}
        }
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
        if(generate_network)
        {
            params.start_day = 1;
        }
        else
        {
            if(pv[19]=="*"){
                std::cout << "ERROR (config 19): No start day option given in config file." << std::endl;}
            params.start_day_option = stringToNum<int>(pv[19]);
            if(params.start_day_option == 0) {params.start_day = rand_int(1, 366);}
            else {params.start_day = params.start_day_option;}
        }

        //What disease to simulate.
        if(pv[20] == "*"){
            std::cout << "ERROR (config 20): Must chose whether to simulate FMD (0) or BTB (1). No other options are allowed." << std::endl; exitflag=1;}
        int inf_type_option = stringToNum<int>(pv[20]);
        switch(inf_type_option)
        {
            case 0:
                params.infectionType = InfectionType::FMD;
                break;
            case 1:
                params.infectionType = InfectionType::BTB;
                break;
            default:
                params.infectionType = InfectionType::FMD;
                std::cout << "ERROR (config 20): Must chose whether to simulate FMD (0) or BTB (1). No other options are allowed." << std::endl; exitflag=1;
                break;
        }

        if(params.infectionType == InfectionType::BTB and !generate_network){
            params.days_per_timestep = 30; }
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

		//Market behavior with regards to WH-transmission.
		if(pv[23]=="*"){
            std::cout << "ERROR (config 23): Market within-herd transmission option must be set." << std::endl;
            exitflag=1;
		}
		if(params.infectionType == InfectionType::BTB)
        {
            int market_WH_transmission = std::stoi(pv[23]);
            switch(market_WH_transmission)
            {
            case 0:
                params.within_spread_at_markets = true; break;
            case 1:
                params.within_spread_at_markets = false; break;
            default:
                std::cout << "ERROR (config 23): Market within-herd transmission option must be 0, 1." << std::endl;
                exitflag=1;
                break;
            }
        }
        else
        {
            params.market_shipment_exposure_prob = std::stod(pv[23]);
        }


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

		// Partial transition off (0), or on (2). Only affects FMD simulations, always on for bTB.
        if(params.infectionType == InfectionType::BTB)
        {
            params.partial = 1;
        }
        else
        {
            if (pv[33]!="0" && pv[33]!="1"){std::cout << "ERROR (config 33): Partial transition flag must be 0 (off) or 1 (on)." << std::endl; exitflag=1;}
            params.partial=stringToNum<int>(pv[33]);
        }

        // Partial transition parameter values
        if (params.partial == 1){
            params.partialParams = stringToNumVec(pv[34]);
            checkExit = checkPositive(params.partialParams, 34); if (checkExit==1){exitflag=1;}
            if (params.partialParams.size() != 5){
                std::cout<<"ERROR (34): Expecting 5 parameter values for partial transition extension." <<std::endl; exitflag=1;
            }
            if(std::get<1>(params.latencyParams)!=0){
                 std::cout<<"ERROR (31): Variance of days from premises exposure to infectiousness must be 0 if partial transistion flag is on." <<std::endl; exitflag=1;
            }
        }

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
                params.shipment_kernel = "USAMMv1"; //USAMMv1
                params.usamm_version = 1;
                break;
            case 2:
                params.shipment_kernel = "USAMMv3"; //USAMMv3
                params.usamm_version = 3;
                break;
            case 3:
                params.shipment_kernel = "power_exp"; //USAMMv2 kernel 1
                params.usamm_version = 2;
                break;
            case 4:
                params.shipment_kernel = "one_minus_exp"; //USAMMv2 kernel 2
                params.usamm_version = 2;
                break;
            case 5:
                params.shipment_kernel = "local"; //USAMMv2 kernel 3
                params.usamm_version = 2;
                break;
            default:
                params.shipment_kernel = "::UNKNOWN, CHECK CONFIG FILE::";
                params.usamm_version = -1;
                break;
            }
            if(params.infectionType == InfectionType::BTB and
               params.shipments_on == 1 and params.usamm_version != 3)
            {
                std::cout << "When simulating bTB infection, shipments must either be set to off or to use USAMMv3 (shipment option = 2)." << std::endl; exitflag = 1;
            }

            if (pv[42]!="*")
            {
                params.shipment_rate_factor = std::stod(pv[42]);
                if(params.shipment_rate_factor < 0.0){ std::cout << "ERROR (config 42): The shipment rate scaling factor must be positive." << std::endl; exitflag=1; }
            }
            //I've removed the option to use different shipment generation methods at different times in the sim as I believe it was unused /S.
//            if (pv[42]=="*"){std::cout << "ERROR (config 42): No shipment method start time(s) specified." << std::endl; exitflag=1;}
//            params.shipMethodTimeStarts = stringToIntVec(pv[42]);
//            if (params.shipMethods.size() != params.shipMethodTimeStarts.size()){
//                std::cout << "ERROR (config 41-42): Number of methods and start times provided must be the same ("<<params.shipMethods.size()<<" method(s) and "<<params.shipMethodTimeStarts.size()<<" start time(s) provided)."<<std::endl; exitflag=1;}
//            if (params.shipMethodTimeStarts.at(0) != 1){std::cout << "ERROR (config 42): First county-level shipment method start time must be 1."<<std::endl; exitflag=1;}
             // check that each value is between 1 and number of timesteps

//            bool rangeflag = 0;
//            for (auto& t:params.shipMethodTimeStarts){if(t<1 || t>params.timesteps){rangeflag=1;}}
//            if (rangeflag){std::cout << "Warning (config 42): Simulation timespan (config 13) is shorter than largest timepoint for county-level shipment methods - some methods may not be used." << std::endl;}
//            params.shipMethodTimeStarts.emplace_back(params.timesteps+1); // add this b/c whichElement function needs a maximum (element won't be used)

            if(params.shipMethods > 0 and pv[44]=="*"){std::cout << "ERROR (config 44): No USAMM parameter file(s) specified." << std::endl; exitflag=1;}
            params.USAMM_parameter_files = stringToStringVec(pv[44]);
            if(pv[44]!="*" and params.USAMM_parameter_files.size() != params.species.size()){std::cout << "ERROR (config 12 & 44): The number of USAMM parameter files and number of species in the config file does not match." << std::endl; exitflag=1;}
            if(pv[44]!="*" and params.shipments_on == 1)
            {
                params.USAMM_post_lengths.reserve(params.USAMM_parameter_files.size());
                for(std::string usamm_post_fn : params.USAMM_parameter_files)
                {
                    std::ifstream usamm_post_f(usamm_post_fn);
                    if(usamm_post_f.is_open())
                    {
                        params.USAMM_post_lengths.push_back(get_n_lines(usamm_post_f));
                        usamm_post_f.close();
                    }
                    else{ std::cout << "ERROR (config 44): Failed to open USAMM posterior distribution file " << usamm_post_fn << std::endl; exitflag=1; }
                }
            }
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
            params.USAMM_temporal_n_timesteps = params.USAMM_temporal_n_days; //For fmd timestep=1 day
            if(params.infectionType == InfectionType::BTB and !generate_network)
            {
                size_t n_periods = params.USAMM_temporal_n_days.size();
                params.USAMM_temporal_n_timesteps.clear();
                params.USAMM_temporal_n_timesteps.resize(n_periods, 12/n_periods); //For BTB there are three months / quarter and one month is a timestep so 3 timesteps per period. This will fail if a period is not something divisible by 12 months, which it shouldn't be.
            }


//            if(params.shipMethods.front() == 2 and pv[47]=="*"){std::cout << "ERROR (config 47): No USAMM origin covariate file(s) specified." << std::endl; exitflag=1;}
            params.USAMM_ocov_files = stringToStringVec(pv[47]);
//            if(params.shipMethods.front() == 2 and pv[48]=="*"){std::cout << "ERROR (config 48): No USAMM destination covariate file(s) specified." << std::endl; exitflag=1;}
            params.USAMM_dcov_files = stringToStringVec(pv[48]);
//            if(params.shipMethods.front() == 2 and pv[49]=="*"){std::cout << "ERROR (config 49): No USAMM supernode files specified." << std::endl; exitflag=1;}
//            params.USAMM_supernode_files = stringToStringVec(pv[49]);
            params.USAMM_supernode_files = std::vector<std::string>(params.species.size(), "*"); //Turned off for now, but left in code if there is interest for this in the future.

            if (params.USAMM_ocov_files.size() != params.species.size()){
                std::cout<<"ERROR (config 12 & 47): Different numbers of origin covariate input files: "<<params.species.size()<<" species and "
                         <<params.USAMM_ocov_files.size()<<" origin covariate files." <<std::endl; exitflag=1;}
            if (params.USAMM_dcov_files.size() != params.species.size()){
                std::cout<<"ERROR (config 12 & 47): Different numbers of destination covariate input files: "<<params.species.size()<<" species and "
                         <<params.USAMM_dcov_files.size()<<" destination covariate files." <<std::endl; exitflag=1;}
            if (params.USAMM_supernode_files.size() != params.species.size()){
                std::cout<<"ERROR (config 12 & 47): Different numbers of supernode input files: "<<params.species.size()<<" species and "
                         <<params.USAMM_supernode_files.size()<<" supernode files." <<std::endl; exitflag=1;}


			params.statuses_to_generate_shipments_from = { "inf" };
			params.statuses_to_generate_slaughter_shipments_from = { "inf" };
			std::string exposed_shipment_option = pv[50];
			if(params.infectionType == InfectionType::FMD and exposed_shipment_option.compare("1") == 0)
            {
                params.statuses_to_generate_shipments_from.push_back("exp");
                params.statuses_to_generate_slaughter_shipments_from.push_back("exp");
            }
            else if(params.infectionType == InfectionType::BTB)
            {
                params.statuses_to_generate_shipments_from = { "btb" };
                params.statuses_to_generate_slaughter_shipments_from = { "btb" };
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
		//Flexible control implimentation
		params.flexibleControl = pv[62];
	//	std::string options_flex = "flex,*";
			if (verbose>1){if(pv[62].compare("flex")==0){std::cout<< "Flexible control option initiated."<<std::endl;
						  }else{std::cout<<"Flexible control option off."<<std::endl;}}
		//	options = stringToStringVec(options_flex);
		//	std::vector<std::string> flexoption = stringToStringVec(pv[62]);
		//	checkExit = limitedValues(flexoption, options, 62); if (checkExit==1){exitflag=1;} //gets in the way...
		// Control rules
		std::vector<std::string> triggers = stringToStringVec(pv[63]);
			std::string options_triggers = "newPremReportsOverX,newRegionReportsOverX";
			options = stringToStringVec(options_triggers);
		//	checkExit = limitedValues(triggers, options, 63); if (checkExit==1){exitflag=1;}
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
			std::string options_orders = "earliest,largest,closest,farthest,smallest"; // can add "largest, closest" later
			options = stringToStringVec(options_orders);
			checkExit = limitedValues(orders, options, 51); if (checkExit==1){exitflag=1;}
			if (orders.size() != triggers.size()){
			std::cout << "ERROR (config 63 and 67): Different number of arguments and control triggers." << std::endl;
			exitflag=1;
			}
		std::string fCFiles = pv[68];
		if (verbose>1){std::cout<< "Flexible control file: " << pv[68]<<std::endl;}

		//pv[69]... pv[70]
			// Put control parameters into proper containers
			for (size_t ct2 = 0; ct2 < params.controlTypes.size(); ++ct2){
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

			
			params.flexibleControlFiles = fCFiles;
			
			//Put flexible control parameters from "flex.txt" into proper controlRule
			if ((params.flexibleControl).compare("flex") == 0){

				if (verbose>1){std::cout << "FM: entering flexibleControlUpdateFile" << std::endl;}

				std::string controlUpdateFile = params.flexibleControlFiles;

					std::ifstream f(controlUpdateFile);
					if(!f){std::cout << "FM: Flexible control file not found. Check config 68. Exiting..." << std::endl; exit(EXIT_FAILURE);}

					if(f.is_open())
					{
					  skipBOM(f);
				    if (verbose>0){std::cout << "FM: Flexible control file open, loading thresholds..." << std::endl;}

						while(! f.eof())
						{
                            std::string line;
							getline(f, line); // get line from file "f", save as "line"
							std::vector<std::string> line_vector = split(line, '\t'); // separated by space
								if(! line_vector.empty()) // if line_vector has something in it
								{
                                    for (size_t cr = 0; cr < (line_vector.size()); cr+=6){  //could generalize this line of code, but is fully functional
                                        if (line_vector[cr]!="percentIncrease" && line_vector[cr+5]!="0") {
                                            std::cout << "ERROR (config 68): Flexible control file has incorrect format. Problematic line is: "<< line_vector[cr]<< " and " << line_vector[cr+5] <<". \n Set decision delay equal to 0 when percentIncrease is not the control trigger. Exiting..." << std::endl;

                                            exit(EXIT_FAILURE);
                                        }
                                        params.controlRules.emplace_back(					//with current flexible control implementation.
											controlRule{									//Works with a flex.txt file with 6 columns: trigger type, threshold, control action,
												line_vector[cr],							//targets, prioritization, decision delay.
												stringToNum<double>(line_vector[cr+1]),
												line_vector[cr+2],
												stringToNum<double>(line_vector[cr+3]),
												line_vector[cr+4],
												stringToNum<double>(line_vector[cr+3])*stringToNum<double>(line_vector[cr+3]),
                                                stringToNum<double>(line_vector[cr+5])
											});

										if(verbose>1){std::cout<<' '<<line_vector[cr];
											std::cout<<' '<<line_vector[cr+1];
											std::cout<<' '<<line_vector[cr+2];
											std::cout<<' '<<line_vector[cr+3];
											std::cout<<' '<<line_vector[cr+4];
                                            std::cout<<' '<<line_vector[cr+5]<<std::endl;
										}

										if (stringToNum<double>(line_vector[cr+3]) == -1){
											params.dangerousContacts_on = 1;
											params.dcControlTypes.emplace_back(line_vector[cr+2]); // Put control rules into containers
										}
									}
								} // end if line has something in it
					} // end while not end of file
				} // if file is open
				if(verbose > 1){std::cout << "Flexible control has " << (params.controlRules).size() << " steps." << std::endl;}
			} else {
				for (size_t cr = 0; cr < triggers.size(); cr++){
				params.controlRules.emplace_back(
					controlRule{
						triggers[cr],
						thresholds[cr],
						actions[cr],
						targets[cr],
						orders[cr],
						targets[cr]*targets[cr], // radius squared, stored pre-calculated for speed
                        0 // decision delay initialized with 0. Only used during flexible control strategy runs.
                    });
					if (targets[cr] == -1){
						params.dcControlTypes.emplace_back(actions[cr]); // Put control rules into containers
					}
				}
			}
		}


// diagnostics

//should infected prems be tested prior to being reported
	params.testSuspects_on = true;
 if (pv[78] == "*"){
	params.testSuspects_on = false;}

//is diagnostics on?
params.diagnostic_on = true;
if (pv[79]=="*" || pv[91]=="*"){
params.diagnostic_on = false;}

if(params.diagnostic_on == true){
		params.diagnosticTypes = stringToStringVec(pv[79]);
			std::string options_dt = "pcr,elisa,penSide,series,diva";
			std::vector<std::string> options = stringToStringVec(options_dt);
			checkExit = limitedValues(params.diagnosticTypes, options, 79); if (checkExit==1){exitflag=1;}
		// Diagnostic - constraint function names
		std::vector<std::string> dFuncs = stringToStringVec(pv[80]);
			if (dFuncs.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 80): Different number of arguments and diagnostic types." << std::endl;
			exitflag=1;
			}
			std::string options_df = "noLimit,diagnosticDailyLimit"; //nationalLimit and stateSum can be added
			options = stringToStringVec(options_df);
			checkExit = limitedValues(dFuncs, options, 80); if (checkExit==1){exitflag=1;}
		// Diagnostic - constraint function parameters
		std::vector<std::string> dFuncParams = semicolonStringToStringVec(pv[81]);
			if (dFuncParams.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 81): Different number of arguments and diagnostic types. Provide 0 for diagnosticConstraintFunctions = 'noLimit'. Or are parameter sets not separated by semicolon?" << std::endl;
			exitflag=1;
			}
		// Diagnostic - spatial scale for each type
		std::vector<std::string> diagnosticSpatialScales = stringToStringVec(pv[82]);
			std::string diagnostic_options_scales = "premises"; //county or state can be added
			options = stringToStringVec(diagnostic_options_scales);
			checkExit = limitedValues(diagnosticSpatialScales, options, 82); if (checkExit==1){exitflag=1;}
			if (diagnosticSpatialScales.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 82): Different number of arguments and diagnostic types." << std::endl;
			exitflag=1;
			}
		// Diagnostic - external files relating to constraint functions
		std::vector<std::string> dFuncFiles = semicolonStringToStringVec(pv[83]);
			if (dFuncFiles.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 83): Different number of arguments and diagnostic types. Are sets not separated by semicolon?" << std::endl;
			exitflag=1;
			}
		// Diagnostic - types for each external file relating to constraint functions
		std::vector<std::string> dFuncFileTypes = semicolonStringToStringVec(pv[84]);
			std::string diagnostic_options_filetypes = "NA"; //currently not set up to take files. Can be added
			options = stringToStringVec(diagnostic_options_filetypes);
			checkExit = limitedValues(dFuncFileTypes, options, 84); if (checkExit==1){exitflag=1;}
			if (dFuncFileTypes.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 84): Different number of arguments and diagnostic types. Are sets not separated by semicolon?" << std::endl;
			exitflag=1;
			}
		// Diagnostic - time lags for each diagnostic test. Lag covers time from implementation to completion
		std::vector<double> testStartToCompleteMeans = stringToNumVec(pv[85]);
			checkExit = checkPositive(testStartToCompleteMeans, 85); if (checkExit==1){exitflag=1;}
			if (testStartToCompleteMeans.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 85): Different number of arguments and diagnostic types." << std::endl;
			exitflag=1;
			}
		// Diagnostic - time lags for each diagnostic test. Lag covers time from implementation to completion
		std::vector<double> testStartToCompleteVars = stringToNumVec(pv[86]);
			checkExit = checkPositive(testStartToCompleteVars, 86); if (checkExit==1){exitflag=1;}
			if (testStartToCompleteVars.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 86): Different number of arguments and diagnostic types." << std::endl;
			exitflag=1;
			}
		// Diagnostic - sensitivity of each diagnostic: mean
		std::vector<double> sensAllMeans = stringToNumVec(pv[87]);
		checkExit = checkPositive(sensAllMeans, 87); if (checkExit==1){exitflag=1;}
		if (sensAllMeans.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 87): Different number of arguments and diagnostic types." << std::endl;
			exitflag=1;
			}
		// Diagnostic - sensitivity of each diagnostic: variance
		std::vector<double> sensAllVars = stringToNumVec(pv[88]);
		checkExit = checkPositive(sensAllVars, 88); if (checkExit==1){exitflag=1;}
		if (sensAllVars.size() != params.diagnosticTypes.size()){
			std::cout << "ERROR (config 79 and 88): Different number of arguments and diagnostic types." << std::endl;
			exitflag=1;
			}
		//pv[89]...pv[90] could add specificity here
		// Diagnostic rules
		std::vector<std::string> diagnostic_triggers = stringToStringVec(pv[91]);
			std::string diagnostic_options_triggers = "newPremReportsOverX,newPremSuspectsOverX,newPremVaxsOverX";
			options = stringToStringVec(diagnostic_options_triggers);
			checkExit = limitedValues(diagnostic_triggers, options, 91); if (checkExit==1){exitflag=1;}

		std::vector<double> diagnostic_thresholds = stringToNumVec(pv[92]);
			checkExit = checkPositive(diagnostic_thresholds, 92); if (checkExit==1){exitflag=1;}
			if (diagnostic_thresholds.size() != diagnostic_triggers.size()){
			std::cout << "ERROR (config 91 and 92): Different number of arguments and diagnostic triggers." << std::endl;
			exitflag=1;
			}
		std::vector<std::string> diagnostic_actions = stringToStringVec(pv[93]);
			options = stringToStringVec(options_dt); // diagnostic types from line 79
			checkExit = limitedValues(diagnostic_actions, options, 93); if (checkExit==1){exitflag=1;}
			if (diagnostic_actions.size() != diagnostic_triggers.size()){
			std::cout << "ERROR (config 91 and 93): Different number of arguments and diagnostic triggers." << std::endl;
			exitflag=1;
			}
		std::vector<double> diagnostic_targets = stringToNumVec(pv[94]);	// should be -1, 0, or positive.
			if (diagnostic_targets.size() != diagnostic_triggers.size()){
			std::cout << "ERROR (config 91 and 94): Different number of arguments and diagnostic triggers." << std::endl;
			exitflag=1;
			}
			for (auto&t:diagnostic_targets){
				if (t==-1){params.dangerousContacts_on = 1;}
			}

// Should also check that any regional diagnostic types (scale != premises) can only have target=0
		std::vector<std::string> diagnostic_orders = stringToStringVec(pv[95]);
			std::string diagnostic_options_orders = "earliest"; // can add "largest, closest" later or other attribute of premises
			options = stringToStringVec(diagnostic_options_orders);
			checkExit = limitedValues(diagnostic_orders, options, 95); if (checkExit==1){exitflag=1;}
			if (diagnostic_orders.size() != diagnostic_triggers.size()){
			std::cout << "ERROR (config 79 and 95): Different number of arguments and diagnostic triggers." << std::endl;
			exitflag=1;
			}

		//pv[95]... pv[100]

		// exposed to suspect lag (only active when one of the targets for diagnostics is 0 or IP)- index
		checkExit = checkMeanVar(pv[100],100,"index suspected"); if (checkExit==1){exitflag=1;}
		tempVec = stringToNumVec(pv[100]);
		params.indexSuspectLag = std::make_tuple(tempVec[0],tempVec[1]);
		// exposed to suspect lag (only active when one of the targets for diagnostics is 0 or IP)- not index
		checkExit = checkMeanVar(pv[101],101,"non-index suspected"); if (checkExit==1){exitflag=1;}
		tempVec = stringToNumVec(pv[101]);
		params.nonDCSuspectLag = std::make_tuple(tempVec[0],tempVec[1]);
		// exposed to suspect lag for DCs exposed after DC status change
		checkExit = checkMeanVar(pv[102],102,"DC suspected"); if (checkExit==1){exitflag=1;}
		tempVec = stringToNumVec(pv[102]);
		params.dcSuspectLag = std::make_tuple(tempVec[0],tempVec[1]);

		//pv[103]... pv[108]

		// Put diagnostic parameters into proper containers
		for (size_t dt2 = 0; dt2 < params.diagnosticTypes.size(); ++dt2){
			std::string d_type = params.diagnosticTypes.at(dt2);

			params.diagnosticScales[d_type] = diagnosticSpatialScales.at(dt2);
			params.diagnosticConstraintFunctions[d_type] = dFuncs.at(dt2);
			params.diagnosticConstraintFuncParams[d_type] = stringToNumVec(dFuncParams.at(dt2));
			params.diagnosticConstraintFuncFiles[d_type] = stringToStringVec(dFuncFiles.at(dt2));
			params.diagnosticConstraintFuncFileTypes[d_type] = stringToStringVec(dFuncFileTypes.at(dt2));
			params.testStartToCompleteLag[d_type] = std::make_tuple(testStartToCompleteMeans.at(dt2),
				testStartToCompleteVars.at(dt2));
			params.sensitivity[d_type] = std::make_tuple(sensAllMeans.at(dt2),
				sensAllVars.at(dt2));
		}	// end "for each diagnostic type"

		// Put diagnostic rules into containers
		for (size_t dr = 0; dr < diagnostic_triggers.size(); dr++){
			params.diagnosticRules.emplace_back(
				diagnosticRule{
					diagnostic_triggers[dr],
					diagnostic_thresholds[dr],
					diagnostic_actions[dr],
					diagnostic_targets[dr],
					diagnostic_orders[dr],
					diagnostic_targets[dr]*diagnostic_targets[dr] // radius squared, stored pre-calculated for speed
				});
			if (diagnostic_targets[dr] == -1){
				params.dcDiagnosticTypes.emplace_back(diagnostic_actions[dr]);
			}
		}
} // end "if diagnostics_on==true"

if(params.control_on == true || params.diagnostic_on ==true){
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
					std::vector<std::string> options = stringToStringVec(options_dc);
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
					params.dcRiskScale[params.StringToDiseaseStatus.at(statusAtSourceReport)] = dcScale;
					// save largest of scale values
					if (dcScale > params.maxDCScale){
						params.maxDCScale = dcScale;
					}
 				} // end "if 2 arguments provided in set"
			} // end "for each set of arguments"
		} // end "if dangerousContacts_on"
	}


    if(params.infectionType == InfectionType::BTB)
    {
        if(params.diagnostic_on ==true)
        {
            if(pv[110] == "*"){
                std::cout << "ERROR (110): Simulating bTB with diagnostics requires a slaughtershed file." << std::endl;
                exitflag = 1;
            } else {
                params.btb_slaughtershed_fname = pv[110];
            }
            if(pv[111] == "*"){
                std::cout << "ERROR (111): Simulating bTB with diagnostics requires a slaughter proportions file." << std::endl;
                exitflag = 1;
            } else {
                params.btb_slaughter_prop_fname = pv[111];
            }
        }
        //BTB infection parameters [120 - 129]
        if(pv[122] == "*") {
            std::cout << "ERROR (122): Simulating bTB requires a set of local spread kernel parameters." << std::endl;
            exitflag = 1;
        } else {
            std::vector<std::string> temp_v = semicolonStringToStringVec(pv[122]); //One vector / quarter.
            for(auto v : temp_v)
            {
                params.btb_local_kernel_parameters.push_back(stringToNumVec(v));
            }

        }

        if(pv[123] == "*") {
            std::cout << "ERROR (123): Simulating bTB requires a set of wildlife transmission kernel parameters." << std::endl;
            exitflag = 1;
        } else {
            std::vector<std::string> temp_v = semicolonStringToStringVec(pv[123]); //One vector / quarter.
            for(auto v : temp_v)
            {
                params.btb_wildlife_kernel_parameters.push_back(stringToNumVec(v));
            }
        }

        if(pv[124] == "*" or stringToNum<double>(pv[124]) < 0.0 or stringToNum<double>(pv[124]) > 1.0) {
            std::cout << "ERROR (124): Simulating bTB requires the relative importance of the wildlife kernel component to be specified as a number between 0.0 and 1.0." << std::endl;
            exitflag = 1;
        } else {
            params.btb_rel_wildlife_kernel_weight = stringToNum<double>(pv[124]);
        }

//        if(pv[125] == "*" or stringToNum<double>(pv[125]) < 0.0) {
//            std::cout << "ERROR (125): Simulating bTB requires a scaling factor of the transmissibility of the mixed local spread kernel to be specified as a >= 0.0." << std::endl;
//            exitflag = 1;
//        } else {
//            params.btb_mixed_kernel_scaling_factor = stringToNum<double>(pv[125]);
//        }

        //Removals of infected animals at the end of timesteps.
//		if(pv[128]=="*"){
//            std::cout << "ERROR (config 128): Option for end of time step removal of infected from markets must be set." << std::endl;
//            exitflag=1;
//		}
//		int market_wipe_option = std::stoi(pv[128]);
//		switch(market_wipe_option)
//		{
//        case 0:
//            params.wipe_markets_each_t = false; break;
//        case 1:
//            params.wipe_markets_each_t = true; break;
//        default:
//            std::cout << "ERROR (config 128): Option for end of time step removal of infected from markets must be 0, 1." << std::endl;
//            exitflag=1;
//            break;
//		}

        if(params.infectionType == InfectionType::BTB and pv[129] == "*") {
            std::cout << "ERROR (129): Simulating bTB requires a wildlife density file to be provided." << std::endl;
            exitflag = 1;
        } else {
            params.btb_wildlife_density_file = pv[129];
        }

//        BTB within-herd model parameters [130-149]
//        if(pv[130] == "*") {
//            std::cout << "ERROR (130): Simulating bTB requires a total yearly export rate to be specified." << std::endl;
//            exitflag = 1;
//        } else {
//
//        }
        params.btbWithinHerdParams.export_rate = 0.0; //Shipments out of infected premises handled by USAMM + slaughter surveillance piece.

//        if(pv[131] == "*") {
//            std::cout << "ERROR (131): Simulating bTB requires a proportion of the total export that takes place in quarter 4 to be specified." << std::endl;
//            exitflag = 1;
//        } else {
//            params.btbWithinHerdParams.export_prop_q4 = stringToNum<double>(pv[131]);
//        }
        params.btbWithinHerdParams.export_prop_q4 = 0.0; //Shipments out of infected premises handled by USAMM + slaughter surveillance piece.

        if(pv[132] == "*") {
            std::cout << "ERROR (132): Simulating bTB requires a birth rate of dairy farms to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.dairy_birth_rate = stringToNum<double>(pv[132]);
        }

        if(pv[133] == "*") {
            std::cout << "ERROR (133): Simulating bTB requires a mortality rate to be specified." << std::endl;
            exitflag = 1;
        } else {
            std::vector<double> tempv = stringToNumVec(pv[133]);
            params.btbWithinHerdParams.beef_mortality_rate = tempv[0];
            params.btbWithinHerdParams.dairy_mortality_rate = tempv[1];
        }

        if(pv[134] == "*") {
            std::cout << "ERROR (134): Simulating bTB requires an import rate to be specified." << std::endl;
            exitflag = 1;
        } else {
            std::vector<double> tempv = stringToNumVec(pv[134]);
            params.btbWithinHerdParams.beef_import_rate = tempv[0];
            params.btbWithinHerdParams.dairy_import_rate = tempv[1];
        }

        if(pv[135] == "*") {
            std::cout << "ERROR (135): Simulating bTB requires a within-herd transmission rate to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.transmission_rate = stringToNum<double>(pv[135]);
        }

        if(pv[136] == "*") {
            std::cout << "ERROR (136): Simulating bTB requires a cattle-to-cattle contact rate to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.cattle_cattle_contact_rate = stringToNum<double>(pv[136]);
        }

//        if(pv[137] == "*") {
//            std::cout << "ERROR (137): Simulating bTB requires a cattle-to-wildlife contact rate to be specified for each quarter (comma separated)." << std::endl;
//            exitflag = 1;
//        } else {
//            std::vector<double> tempv = stringToNumVec(pv[137]);
//            params.btbWithinHerdParams.cattle_wildl_contact_rate_q1 = tempv[0];
//            params.btbWithinHerdParams.cattle_wildl_contact_rate_q2 = tempv[1];
//            params.btbWithinHerdParams.cattle_wildl_contact_rate_q3 = tempv[2];
//            params.btbWithinHerdParams.cattle_wildl_contact_rate_q4 = tempv[3];
//        }
        //Removed since it's an external/wildlife source and thats already covered.
        params.btbWithinHerdParams.cattle_wildl_contact_rate_q1 = 0.0;
        params.btbWithinHerdParams.cattle_wildl_contact_rate_q2 = 0.0;
        params.btbWithinHerdParams.cattle_wildl_contact_rate_q3 = 0.0;
        params.btbWithinHerdParams.cattle_wildl_contact_rate_q4 = 0.0;

        if(pv[138] == "*") {
            std::cout << "ERROR (138): Simulating bTB requires a mean and rate of gamma distributed transition rate from class E1 to E2U or E2R to be specified (comma-separated)." << std::endl;
            exitflag = 1;
        } else {
            std::vector<double> tempv = stringToNumVec(pv[138]);
            params.btbWithinHerdParams.sigma_1_mean = tempv[0];
            params.btbWithinHerdParams.sigma_1_rate = tempv[1];
        }

        if(pv[139] == "*") {
            std::cout << "ERROR (139): Simulating bTB requires a meean and rate of gamma distributed transition rate from class E2R to IR to be specified (comma-separated)." << std::endl;
            exitflag = 1;
        } else {
            std::vector<double> tempv = stringToNumVec(pv[139]);
            params.btbWithinHerdParams.sigma_2_mean = tempv[0];
            params.btbWithinHerdParams.sigma_2_rate = tempv[1];
        }

        if(pv[140] == "*") {
            std::cout << "ERROR (140): Simulating bTB requires an adjustment factor of how many times longer it takes to transition from E2U compared to the time it takes to transition from E2R to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.phi_sigma_2 = stringToNum<double>(pv[140]);
        }

        if(pv[141] == "*") {
            std::cout << "ERROR (141): Simulating bTB requires an adjustment factor of how many times longer it takes for an animal in E2U to become reactive compared to the time it takes to transition from E1 to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.phi_delta_1 = stringToNum<double>(pv[141]);
        }

        if(pv[142] == "*") {
            std::cout << "ERROR (142): Simulating bTB requires an adjustment factor of how many times longer it takes for an animal in IU to become reactive compared to the time it takes to transition from E2R to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.phi_delta_2 = stringToNum<double>(pv[142]);
        }

        if(pv[143] == "*") {
            std::cout << "ERROR (143): Simulating bTB requires an adjustment factor of the probability of sucessfully detecting an exposed animal that mounts an immune response (leaves E1U) when testing (goes into E2R) to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.response_prob_1 = stringToNum<double>(pv[143]);
        }

        if(pv[144] == "*") {
            std::cout << "ERROR (144): Simulating bTB requires an adjustment factor of the probability of sucessfully detecting an unreactive exposed animal with immune response (E2U), meaning it moves into E2R to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.response_prob_2 = stringToNum<double>(pv[144]);
        }

        if(pv[145] == "*") {
            std::cout << "ERROR (145): Simulating bTB requires an adjustment factor of the probability of sucessfully detecting an infectious animal that is not reactive to testing (IU, goes into IR) to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.response_prob_3 = stringToNum<double>(pv[145]);
        }

        if(pv[146] == "*") {
            std::cout << "ERROR (146): Simulating bTB requires a test sensitivity to be specified." << std::endl;
            exitflag = 1;
        } else {
            params.btbWithinHerdParams.test_sensitivity = stringToNum<double>(pv[146]);
        }

        params.btbWithinHerdParams.test_specificity = 0.0;
//        if(pv[147] == "*") {
//            std::cout << "ERROR (147): Simulating bTB requires a test specificity to be specified." << std::endl;
//            exitflag = 1;
//        } else {
//            params.btbWithinHerdParams.test_specificity = stringToNum<double>(pv[147]);
//        }

        params.btbWithinHerdParams.specificity_scaling_factor = 0.0;
//        if(pv[148] == "*") {
//            std::cout << "ERROR (148): Simulating bTB requires a test specificity scaling factor to be specified." << std::endl;
//            exitflag = 1;
//        } else {
//            params.btbWithinHerdParams.specificity_scaling_factor = stringToNum<double>(pv[148]);
//        }

        //Must be called here.
        params.btbWithinHerdParams.update_latency_parameters();
    }


//All checks of parameter vector complete here.
if (exitflag){
	std::cout << "Exiting..." << std::endl;
	exit(EXIT_FAILURE);
}
		// Construct kernel object
		if(params.infectionType == InfectionType::BTB)
        {
            // bTB kernel parameters by idx are 0: phi (local), 1: alpha (local), 2: mean (wildlife), 3: stdev (wildlife), 4: x (relative weight of wildlife [0,1]), 5: scaling factor for entire combined kernel applied after normalization.
            //alpha and beta needs to include a factor to translate from km to m.
            std::vector<double> btb_k_q1_pars = { params.btb_local_kernel_parameters.at(0).at(0), params.btb_local_kernel_parameters.at(0).at(1),
                                                  params.btb_wildlife_kernel_parameters.at(0).at(0), params.btb_wildlife_kernel_parameters.at(0).at(1),
                                                  params.btb_rel_wildlife_kernel_weight, params.btb_mixed_kernel_scaling_factor };

            std::vector<double> btb_k_q2_pars = { params.btb_local_kernel_parameters.at(1).at(0), params.btb_local_kernel_parameters.at(1).at(1),
                                                  params.btb_wildlife_kernel_parameters.at(1).at(0), params.btb_wildlife_kernel_parameters.at(1).at(1),
                                                  params.btb_rel_wildlife_kernel_weight, params.btb_mixed_kernel_scaling_factor };

            std::vector<double> btb_k_q3_pars = { params.btb_local_kernel_parameters.at(2).at(0), params.btb_local_kernel_parameters.at(2).at(1),
                                                  params.btb_wildlife_kernel_parameters.at(2).at(0), params.btb_wildlife_kernel_parameters.at(2).at(1),
                                                  params.btb_rel_wildlife_kernel_weight, params.btb_mixed_kernel_scaling_factor };

            std::vector<double> btb_k_q4_pars = { params.btb_local_kernel_parameters.at(3).at(0), params.btb_local_kernel_parameters.at(3).at(1),
                                                  params.btb_wildlife_kernel_parameters.at(3).at(0), params.btb_wildlife_kernel_parameters.at(3).at(1),
                                                  params.btb_rel_wildlife_kernel_weight, params.btb_mixed_kernel_scaling_factor };

            std::vector<std::vector<double>> quarterly_kernel_parametes( { btb_k_q1_pars, btb_k_q2_pars, btb_k_q3_pars, btb_k_q4_pars } );
            params.kernel = new Local_spread(3, quarterly_kernel_parametes);
        }
		else
        {
            switch (params.kernelType)
            {
                case 0:{
                    params.kernel = new Local_spread(0, { params.kernelParams });
                    break;
                }
                case 1:{
                    params.kernel = new Local_spread(1, params.dataKernelFile);
                    break;
                }
                case 2:{
                    params.kernel = new Local_spread(2, { params.kernelParams });
                    break;
                }
            }
        }


if (verbose>0){std::cout<<"Parameter loading complete."<<std::endl;}
if (verbose>0){std::cout<<"Parameter params.dangerousContacts_on is "<<params.dangerousContacts_on<<std::endl;}
if (verbose>0){std::cout<<"Parameter params.testSuspects_on is "<<params.testSuspects_on<<std::endl;}

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




