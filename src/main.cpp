// main.cpp - controls timesteps, initiating various managers, and output
#include <iostream>
#include <ctime>
#include <stdlib.h>

#include "File_manager.h"
#include "Grid_manager.h"
#include "Grid_checker.h"
#include "Shipment_manager.h"
#include "County.h"


int verboseLevel; // global variable determining console output

int main(int argc, char* argv[])
{
    std::clock_t process_start = std::clock();
	std::string cfile; // Config file
	bool gen_shipment_network = false;
	int n_networks = 1;
// Check for command line arguments
	if(argc == 2){cfile = argv[1];}
	else if(argc < 2){
		std::cout << "ERROR: No arguments specified. Either specify config file "
                  << "as single argument to simulate or config --make_network "
                  << "to only generate shipments. " << std::endl
                  << "Exiting..." << std::endl;
		exit(EXIT_FAILURE);
	}
	else if(argc == 3 and std::string("--make_network").compare(argv[2]) == 0)
    {
        cfile = argv[1];
        gen_shipment_network = true;
        std::cout << "Generating 1 shipment network." << std::endl;
    }
    else if(argc == 4 and std::string("--make_network").compare(argv[2]) == 0)
    {
        cfile = argv[1];
        gen_shipment_network = true;
        try
        {
            n_networks = std::stoi(argv[3]);
        }
        catch(std::exception& e)
        {
            std::cout << "Could not in a meaningful way interpret " << argv[3] <<
                         " as the number of networks to generate. " << std::endl <<
                         "Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Generating " << n_networks << " shipment networks." << std::endl;
    }
	else{
		std::cout << "ERROR: Invalid arguments specified. Should only be "
                  << "config file or config file and --make_network to "
                  << "only generate a shipment network." << std::endl
                  <<  "Exiting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	File_manager fm; // construct file_manager object
	fm.readConfig(cfile, gen_shipment_network); // reads config file, creates parameters object, and checks for errors
	const Parameters* p = fm.getParams();

	// make string: batch_date_time
	std::string batchDateTime = p->batchDateTime;
	// get/format current time


	// write parameters from config file to settings_batchname
	std::string settingsOutFile = "runlog.txt";
	// columns are batchDateTime, config lines 1-70, tab-separated
	std::string printString = fm.getSettings(batchDateTime);
	printLine(settingsOutFile,printString);

	// set values for global,
	verboseLevel = p->verboseLevel;

	int verbose = verboseLevel; // override global value for main here if desired
	int timesteps = p->timesteps;

	// Read in farms, determine xylimits
	std::clock_t loading_start = std::clock();
	Grid_manager G(p);
	// get pointers to full list of farms & cells
//	auto allPrems = G.get_allFarms();
	const auto fipsvec = G.get_allCountiesVec();
	const auto allCells = G.get_allCells(); // will be filled when grid is initiated, for now pointer just exists
 	const auto fipsSpeciesMap = G.get_fipsSpeciesMap(); // for shipments
	std::clock_t loading_end = std::clock();

	if(verbose>0){
		std::cout << std::endl << "CPU time for loading premises: "
		<< 1000.0 * (loading_end - loading_start) / CLOCKS_PER_SEC
		<< "ms." << std::endl;
	}

	// Initiate grid...
	if(!gen_shipment_network)
    {
        std::clock_t grid_start = std::clock();
        // if cell file provided, use that
        if (p->cellFile!="*"){
            std::string cellFile = p->cellFile;
            G.initiateGrid(cellFile);} // reading in 730 cells takes ~45 sec
        // else use uniform params
        else if (p->uniformSide>0){
            G.initiateGrid(p->uniformSide);}
        // else use density params
        else {
            G.initiateGrid(p->densityParams.at(0), // max prems per cell
                           p->densityParams.at(1)); // min cell side}
        }

        std::clock_t grid_end = std::clock();
        double gridGenTimeMS = 1000.0 * (grid_end - grid_start) / CLOCKS_PER_SEC;
if(verbose>0){
    std::cout << "CPU time for generating grid: " << gridGenTimeMS << "ms." << std::endl;
}

    } //End make grid

    if(!gen_shipment_network) //Not making a shipment network, running disease simulation.
    {
        Control_manager Control(p, &G); // pass parameters and Grid_manager pointer
        Diagnostic_manager Diagnostic(p, &G); // pass parameters and Grid_manager pointer
        // Get initially infected (seed) premises
        std::vector<std::vector<Farm*>> seedFarmsByRun;
        G.get_seedPremises(seedFarmsByRun); // saves to seedFarmsByRun
if(verbose>0){
	std::cout << seedFarmsByRun.size() << " seed premises generated." << std::endl;
}

	if (seedFarmsByRun.size() == 0 ){
	  std::cout << "ERROR: No valid seed farms provided/located. Exiting... ";
	  exit(EXIT_FAILURE);
	}

        if(verbose>0){
            if(p->infectionType == InfectionType::FMD and p->partial==1) std::cout << "FMD like within herd dynamics will be implemented." << std::endl;
        }

    //~~~~~~~~~~~~~~~~~~ Loop starts here
    for (size_t r=1; r<=seedFarmsByRun.size(); r++){
        std::clock_t rep_start = std::clock();
if(verbose>0){
        std::cout << "Replicate "<< r << " starting." << std::endl;
}
        // load initially infected farms and instantiate Status manager
        // note that initial farms are started as exposed
    	std::vector<Farm*> seedFarms = seedFarmsByRun[r-1];
    	int rep_start_day = p->start_day_option;
        if(rep_start_day == 0) //If option = 0 we set a random start day, otherwise we keep whats in the config file.
        {
            rep_start_day = rand_int(1, 366); //Set a new start day each replicate.
        }
        int day_of_year = rep_start_day;
        int month_of_year = get_month_of_year(day_of_year);
        int current_quarter = get_quarter(day_of_year);
        G.updateQuarterlyFarmSizes(current_quarter);

        Status_manager Status(seedFarms, p, &G, &Control, &Diagnostic); // seeds initial exposures, provides parameters, grid, control, diagnostic
        Shipment_manager Ship(fipsvec, fipsSpeciesMap, &Status, p->species, p); // modify to pass grid manager, p
        Grid_checker gridCheck(allCells, &Status, p);
        // control resources

        int t=0;
        if(p->shipments_on and (p->usamm_version == 1 or p->usamm_version == 2))
        {
            G.initShippingParametersUSAMMv2(1, rep_start_day);
        }
        else if(p->shipments_on and p->usamm_version == 3)
        {
            G.initShippingParametersUSAMMv3(1, rep_start_day);
        }
        std::vector<Farm*> focalFarms; //Stores infectious farms for the local spread component.
        std::vector<Farm*> focalFarmsShipments; //Stores both infectious and exposed premises for the shipment component.
        bool potentialTx = 1;

      int t_in_days = 1; //timestep as days, will be same as t for fmd.
      int last_quarter = 0; // 0 = Jan-Mar, 1 = Apr-Jun, 2=Jul-Sep, 3=Oct-Dec.
      std::vector<Shipment*> rep_shipments; //Just to keep track of all shipments created during the replicate so they can be accessed and deleted.
      rep_shipments.reserve(100000);
      while (t<timesteps && potentialTx){ // timesteps, stop early if dies out
            std::clock_t timestep_start = std::clock();

            ++t; // starts at 1, ends at timesteps
            day_of_year = get_day_of_year(t_in_days, rep_start_day);
            current_quarter = get_quarter(day_of_year); //Returns INDEX of quarter i.e. Jan-Mar=0, Apr-Jun=1, Jul-Sep=2, Oct-Dec=3

            if(t == 1)
            {
                last_quarter = current_quarter;
                G.updateQuarterlyFarmSizes(current_quarter);
            }
if(verbose>1){
            std::cout << std::endl <<"Beginning timestep "<< t
                      << ", " << t_in_days << " days have passed, day of the year is "
                      << day_of_year << ", the current quarter is " << current_quarter
                      << std::endl;
}
            if(p->shipments_on and (p->usamm_version == 1 or
                                    p->usamm_version == 2 or
                                    p->usamm_version == 3))
            {
                G.updateShippingParameters(t, day_of_year);
            }
if(verbose>1){
                std::cout << "Date and temporal shipping parameters updated." << std::endl;
}

            if(p->infectionType == InfectionType::FMD)
            {
                Status.updateDisease(t); // disease status updates
if(verbose>1){
                std::cout << "Disease statuses updated" << std::endl;
}

                Status.updateFile(t); //  file status updates
if(verbose>1){
                std::cout << "File statuses updated" << std::endl;
}
                if (p->diagnostic_on == true){
                    Status.updateDiagnostics(t); // diagnostic status updates
					if(verbose>1){
                		std::cout << "Diagnostic statuses updated" << std::endl;
					}
                }
                if (p->control_on == true){
                    Status.updateControl(t); // control status updates
if(verbose>1){
                    std::cout << "Control statuses updated" << std::endl;
}
                }


                Status.get_premsWithStatus("inf", focalFarms);	// set focalFarms as all farms with disease status inf
                if(current_quarter != last_quarter)
                {
                    //Quarter has changed, so update quarterly farm sizes, for FMD this is just a simple change of premsize.
                    G.updateQuarterlyFarmSizes(current_quarter);
                }
            }
            else if(p->infectionType == InfectionType::BTB)
            {
                if(current_quarter != last_quarter)
                {
                    //Quarter has changed, so update quarterly farm sizes. For BTB this means either add the difference to
                    //the number of susceptible animals or remove the difference equally from all classes and
                    //then update the prevalence on the prem. and if simulating btb,
                    //update the prevalence as well.
                    G.updateQuarterlyFarmSizes(current_quarter);
                    Status.btb_updatePremSizeChangesQuarter(t);
                }
                Status.btb_updateDisease(t, current_quarter); //Move all newly exposed from sus to btb status
                Status.updateFile(t); //  file status updates. Unnecessary for btb?
                if (p->diagnostic_on == true)
                {
                    Status.updateDiagnostics(t); // diagnostic status updates
                }
                if (p->control_on == true)
                {
                    Status.updateControl(t); // control status updates
                }
                Status.get_premsWithStatus("btb", focalFarms); // For btb the focal farms are all with a btb status.
            }

            //Build the set of farms to generate shipments from.
            Status.get_premsWithStatus(p->statuses_to_generate_shipments_from, focalFarmsShipments);


if(verbose>0){
                std::cout <<"Timestep "<<t<<": "
                <<Status.numPremsWithStatus("sus")<<" susceptible, "
                <<Status.numPremsWithStatus("exp")<<" exposed, "
                <<focalFarms.size()<<" infectious, "
                <<Status.numPremsWithStatus("imm")<<" immune premises. "<<std::endl;
                if (p->control_on == true || p->diagnostic_on == true){
                    std::cout << Status.numPremsWithFileStatus("reported") << " reported premises in "
                    <<Status.get_numCountiesReported()<<" counties and "
                    <<Status.get_numStatesReported()<<" state(s)."<<std::endl;
                }
}

                // determine infections that will happen from local diffusion
if(verbose>0){
                std::cout << "Starting grid check (local spread): "<<std::endl;}
                std::clock_t gridcheck_start = std::clock();
                std::vector<Farm*> notSus;
                if(p->infectionType == InfectionType::FMD)
                {
                    // gets newly not-susceptible farms to remove from consideration.
                    // For btb we allow exposed premises to be exposed multiple times from local spread, since we track the animals explicitly.
                    Status.newNotSus(notSus);
                }
                gridCheck.stepThroughCells(focalFarms,notSus, t, current_quarter); // records new exposures via Status_manager::addPremForEval. Same function for FMD and BTB.
                std::clock_t gridcheck_end = std::clock();
                double gridCheckTimeMS = 1000.0 * (gridcheck_end - gridcheck_start) / CLOCKS_PER_SEC;
if(verbose>0){
                std::cout << "CPU time for checking grid: " << gridCheckTimeMS << "ms." << std::endl;
}

                // determine shipments
                std::vector<Shipment*> fs; // fs = farm shipments, where new shipments are saved
                fs.reserve(10000);
                double shipTimeMS = 0.0;
                if(!focalFarmsShipments.empty() and p->shipments_on)
                {
                    std::clock_t ship_start = std::clock();
                    std::string time_period = G.get_time_period(); //Current time period we are in given timestep (i.e Q1. Q2, ...)
                    int usamm_temporal_index = G.get_temporal_index();
                    if(p->usamm_version == 1 or p->usamm_version == 2) //Only generate shipments if there are any farms to generate from.
                    {
                        //This will never happen if simulating btb as that forces the use of USAMMv3
                        size_t days_rem = G.get_rem_days_of_period(); //Number of days remaining of this time period.
                        Ship.makeShipmentsUSAMMv2(t, days_rem, time_period, fs,
                                                  focalFarmsShipments, G.get_farm_types());
                    }
                    else if(p->usamm_version == 3)
                    {
                        size_t day_of_year = get_day_of_year(t, p->start_day);
                        Ship.makeShipmentsUSAMMv3(t, day_of_year, time_period, usamm_temporal_index, fs,
                                                  focalFarmsShipments, G.get_usammv3_parameters_vec());

                    }
                    std::clock_t ship_end = std::clock();
                    shipTimeMS = 1000.0 * (ship_end - ship_start) / CLOCKS_PER_SEC;
                    if(verbose>0){std::cout << "CPU time for shipping: " << shipTimeMS << "ms." << std::endl;}
                    // determine which shipments escape ban and are to susceptible farms
                    Status.filter_shipments(fs, t); // Can output shipping from complete "fs" after this point
                }

                // evaluate if premises-level control prevents any exposures, expose accordingly
                if(p->infectionType == InfectionType::FMD)
                {
                    Status.eval_exposure(t, current_quarter);
                }
                else if(p->infectionType == InfectionType::BTB)
                {
                    Status.eval_slaughterSurveillance(t, current_quarter, r, focalFarmsShipments);
                    Status.btb_evalExposure(t);

                }

                if(verbose>1){std::cout << "SM: eval_exposure complete" << std::endl;}
                Status.add_diagnosticWaitlistMembers(t);
                if(verbose>1){std::cout << "SM: add_diagnosticWaitlistMembers complete" << std::endl;}
                Status.add_started(t);
                if(verbose>1){std::cout << "SM: add_started" << std::endl;}
				Status.eval_testResult(t, current_quarter);
				if(verbose>1){std::cout << "SM: eval_testResult complete" << std::endl;}
                Status.add_waitlistMembers(t);
                if(verbose>1){std::cout << "SM: add_waitlistMembers complete" << std::endl;}
				Status.update_ControlResources(t);
                if(verbose>1){std::cout << "SM: update_ControlResources" << std::endl;}
                Status.add_implemented(t);
                if(verbose>1){std::cout << "SM: add_implemented" << std::endl;}

                // at the end of this transmission day, statuses are now...
                int numSuscept = 0;
                int numExposed = 0;
                size_t numPotentialTx = 0;
                if(p->infectionType == InfectionType::FMD)
                {
                    Status.get_premsWithStatus("inf", focalFarms); // assign "inf" farms as focalFarms
                    numSuscept = Status.numPremsWithStatus("sus");
                    numExposed = Status.numPremsWithStatus("exp");
                    numPotentialTx = focalFarms.size() + numExposed;
                }
                else if(p->infectionType == InfectionType::BTB)
                {
                    Status.btb_updateRecords();
                    numPotentialTx = Status.numPremsWithStatus("btb");
                    numSuscept = Status.numPremsWithStatus("sus") + Status.numPremsWithStatus("rec");

                    //No more transmission this time step. If markets should be wiped, it happens here.
                    //This needs something smarter. With discrete timesteps, this means that markets
                    //doesn't have an effect since any infected animals received this time step will
                    //be removed before they can start to transmit in the next timestep.
//                    if(p->wipe_markets_each_t)
//                    {
//                        Status.wipe_markets();
//                    }
                }

                // write output for details of exposures from this rep, t
                if (p->printDetail > 0){
                    // output detail to file
                    // rep, ID, time, sourceID, route, prevented
                    std::string detOutFile = batchDateTime;
                    detOutFile += "_detail.txt";
                    if (r==1 && t==1){
                        std::string header = "Rep\tExposedID\tatTime\tSourceID\tInfRoute\tControlPrevented\tExposedCounty\tSourceCounty\n";
                        printLine(detOutFile,header);
                    }
                    std::string printString = Status.formatDetails(r,t);
                    printLine(detOutFile, printString);
                }

                //Are there any prems that can still transmit? If not, end FMD sim.
                potentialTx = (numPotentialTx>0 && numSuscept>0);
                if (p->useMaxPrems==1){
                    int totalInf = Status.get_totalPremsWithStatus("inf");
                    if (totalInf > p->maxInfectiousPrems){
                        potentialTx = 0;
                    }
                }

                rep_shipments.insert(rep_shipments.begin(), fs.begin(), fs.end());

                if(p->infectionType == InfectionType::FMD)
                {
                    t_in_days += 1;
                }
                else if(p->infectionType == InfectionType::BTB)
                {
                    t_in_days += p->DAYS_PER_MONTH[month_of_year];
                }
                ++month_of_year;
                if(month_of_year == 12)
                {
                    month_of_year = 0; //Start over with January (0)
                }
                last_quarter = current_quarter;

                std::clock_t timestep_end = std::clock();
                double timestepTimeMS = 1000.0 * (timestep_end - timestep_start) / CLOCKS_PER_SEC;
if(verbose>0){
                std::cout << "CPU time for timestep "<< timestepTimeMS << "ms." << std::endl;
}
        }  	// end "while under time and exposed/infectious and susceptible farms remain"

            std::clock_t rep_end = std::clock();
            double repTimeMS = 1000.0 * (rep_end - rep_start) / CLOCKS_PER_SEC;
            std::cout << "CPU time for batch "<<batchDateTime<<", seed source #"<<r<<" of "
                      <<seedFarmsByRun.size()<<" ("<<t<<" timesteps): " << repTimeMS << "ms." << std::endl;

            if (p->printSummary > 0)
            {
                // output summary to file (rep, days inf, run time)
                // rep, # farms infected, # days of infection, seed farm and county, run time
                std::string sumOutFile = batchDateTime;
                sumOutFile += "_summary.txt";
                if (r==1)
                {
                    std::stringstream ss;
//                    ss << "#\t" << "local_phi.122.a\t" << p->btb_local_kernel_parameters[0][0] << std::endl
//                       << "#\t" << "local_alpha.122.b\t" << p->btb_local_kernel_parameters[0][1] << std::endl
//                       << "#\t" << "rel_wildl.124\t" << p->btb_rel_wildlife_kernel_weight << std::endl
//                       << "#\t" << "dairy_brate.132\t" << p->btbWithinHerdParams.dairy_birth_rate << std::endl
//                       << "#\t" << "b_mortality.133.a\t" << p->btbWithinHerdParams.beef_mortality_rate << std::endl
//                       << "#\t" << "d_mortality.133.b\t" << p->btbWithinHerdParams.dairy_mortality_rate << std::endl
//                       << "#\t" << "b_import.134.a\t" << p->btbWithinHerdParams.beef_import_rate << std::endl
//                       << "#\t" << "d_import.134.b\t" << p->btbWithinHerdParams.dairy_import_rate << std::endl
//                       << "#\t" << "trans_rate.135\t" << p->btbWithinHerdParams.transmission_rate << std::endl
//                       << "#\t" << "cont_rate.136\t" << p->btbWithinHerdParams.cattle_cattle_contact_rate << std::endl
//                       << "#\t" << "lat_adj1.140\t" << p->btbWithinHerdParams.phi_sigma_2 << std::endl
//                       << "#\t" << "lat_adj2.141\t" << p->btbWithinHerdParams.phi_delta_1 << std::endl
//                       << "#\t" << "lat_adj3.142\t" << p->btbWithinHerdParams.phi_delta_2 << std::endl;

                    std::string header = ss.str() + "Rep\tNum_Inf\tnAffCounties\tDuration\tSeed_Farms\tSeed_FIPS\tRunTimeSec\tNum_Reports";
                    if (p->control_on == 1)
                    {
                        for(auto& ct:(p->controlTypes))
                        {
                            std::string ctString = "\t"+ct+"Implemented";
                            ctString += "\t"+ct+"Effective";
                            header+=ctString;
                            for(auto& dct:(p->dcControlTypes))
                            {
                                if (dct==ct)
                                {
                                    std::string dcImp = "\t"+ct+"ImplementedDCSubset";
                                    header+=dcImp;
                                }
                            }

                        }
                    }
                    if (p->diagnostic_on == 1)
                    {
                        for(auto& dt:(p->diagnosticTypes))
                        {
                            std::string dtString = "\t"+dt+"Started";
                            dtString += "\t"+dt+"Complete";
                            header+=dtString;
                            for(auto& dct:(p->dcDiagnosticTypes))
                            {
                                if (dct==dt)
                                {
                                    std::string dcStart = "\t"+dt+"StartedDCSubset";
                                    header+=dcStart;
                                }
                            }

                        }
                    }
                    if (p->dangerousContacts_on == 1)
                    {
                        std::string meandc = "\tmeanDCsPerRP";
                        header+=meandc;
                    }
                    if(p->infectionType == InfectionType::BTB)
                    {
                        std::string n_inf_animals_str = "\tNum_Anim_Inf";
                        std::string n_rec_str = "\tNum_Rec";
                        header += n_inf_animals_str + n_rec_str;
                    }
                    header+="\n";
                    printLine(sumOutFile,header);
                }
                std::string repOut = Status.formatRepSummary(r,t,repTimeMS);
                printLine(sumOutFile,repOut);
            }

            //Delete all the generated shipments.
            for(size_t shp_idx=0; shp_idx<rep_shipments.size(); ++shp_idx)
            {
                delete rep_shipments[shp_idx];
            }

if(verbose>0){
            std::cout << "Replicate "<< r << " complete." << std::endl<< std::endl;
}
    } // end for each replicate loop
} // End if !gen_shipment_network
else if(gen_shipment_network)
{
    if(!p->shipments_on)
    {
        std::cout << "Shipment networks can only be generated if shipments are turned on (config 41)."
                  << " Activate and rerun. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
    for(int netw_number=1; netw_number<=n_networks; netw_number++)
    {
        if(p->usamm_version == 1 or p->usamm_version == 2)
        {
            G.initShippingParametersUSAMMv2(1, p->start_day);
            Shipment_manager Ship(fipsvec, fipsSpeciesMap, p->species, p);
            std::vector<std::string> network_out_files;
            std::string netw_number_s = std::to_string(netw_number);
            for(std::string s : p->species)
            {
                network_out_files.push_back(p->batch + "_" + s + "_" + netw_number_s);
            }
            Ship.makeNetworkUSAMMv2(network_out_files, G);
        }
        else if(p->usamm_version == 3)
        {
            //Make USAMMv3 network here.
            G.initShippingParametersUSAMMv3(1, p->start_day);
            Shipment_manager Ship(fipsvec, fipsSpeciesMap, p->species, p);
            std::vector<std::string> network_out_files;
            std::string netw_number_s = std::to_string(netw_number);
            for(std::string s : p->species)
            {
                network_out_files.push_back(p->batch + "_" + s + "_" + netw_number_s);
            }
            Ship.makeNetworkUSAMMv3(network_out_files, G);
        }
    }
}

    std::clock_t process_end = std::clock();
    double process_time = 1000.0 * (process_end - process_start) / CLOCKS_PER_SEC;
    std::cout << "Entire process was alive for " << process_time << "ms." << std::endl;
    return 0;
} // end main()
