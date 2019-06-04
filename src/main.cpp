// main.cpp - controls timesteps, initiating various managers, and output
#include <iostream>
#include <ctime>
#include <stdlib.h>

#include "File_manager.h"
#include "Grid_manager.h"
#include "Grid_checker.h"
#include "Shipment_manager.h"
#include "County.h"
#include "Population_manager.h"

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
	fm.readConfig(cfile); // reads config file, creates parameters object, and checks for errors
	const Parameters* p = fm.getParams();

	// make string: batch_date_time
	std::string batchDateTime = p->batch;
	// get/format current time
	time_t rawtime;
	struct tm* timeinfo;
	char buffer[80];
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(buffer,80,"_%Y%b%d_%H%M",timeinfo);
	std::string str(buffer);

	batchDateTime += str;
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
	const auto fipsmap = G.get_allCounties();
	const auto allCells = G.get_allCells(); // will be filled when grid is initiated, for now pointer just exists
 	const auto fipsSpeciesMap = G.get_fipsSpeciesMap(); // for shipments
	std::clock_t loading_end = std::clock();

	if(verbose>0){
		std::cout << std::endl << "CPU time for loading premises: "
		<< 1000.0 * (loading_end - loading_start) / CLOCKS_PER_SEC
		<< "ms." << std::endl;
	}

	// Initiate grid...
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

    if(!gen_shipment_network) //Not making a shipment network, running disease simulation.
    {		
        Control_manager Control(p, &G); // pass parameters and Grid_manager pointer
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
            if(p->partial==1) std::cout << "FMD like within herd dynamics will be implemented." << std::endl;
            if(p->partial==2) std::cout << "bTB like within herd dynamics will be implemented." << std::endl;
        }

    //~~~~~~~~~~~~~~~~~~ Loop starts here
    for (int r=1; r<=seedFarmsByRun.size(); r++){
        std::clock_t rep_start = std::clock();
        // load initially infected farms and instantiate Status manager
        // note that initial farms are started as exposed
    	std::vector<Farm*> seedFarms = seedFarmsByRun[r-1];
    	int rep_start_day = p->start_day_option;
        if(rep_start_day == 0) //If option = 0 we set a random start day, otherwise we keep whats in the config file.
        {
            rep_start_day = rand_int(1, 365); //Set a new start day each replicate.
        }


        Status_manager Status(seedFarms, p, &G, &Control); // seeds initial exposures, provides parameters, grid, control
        Shipment_manager Ship(fipsmap, fipsSpeciesMap, &Status, p->shipPremAssignment, p->species, p); // modify to pass grid manager, p
        Grid_checker gridCheck(allCells, &Status, p);
        // control resources

        Population_manager Pop(&Status, p);//initialise Pop manager

        if(p->partial==2){ //if btb infection
            for (auto& f:seedFarms){
                Pop.set_initialFarmSize(f); //set the initial farm size for the seed farms
            }
        }

        int t=0;
        if(p->shipments_on)
        {
            G.initShippingParameters(1, rep_start_day);
        }
        std::vector<Farm*> focalFarms; //Stores infectious farms for the local spread component.
        std::vector<Farm*> focalFarmsShipments; //Stores both infectious and exposed premises for the shipment component.
        bool potentialTx = 1;

      while (t<timesteps && potentialTx){ // timesteps, stop early if dies out
            std::clock_t timestep_start = std::clock();

            ++t; // starts at 1, ends at timesteps
if(verbose>1){
            std::cout << std::endl <<"Beginning timestep "<< t << std::endl;
}
            if(p->shipments_on)
            {
                G.updateShippingParameters(t);
            }
if(verbose>1){
                std::cout << "Date and temporal shipping parameters updated." << std::endl;
}
                Status.updateDisease(t); // disease status updates
if(verbose>1){
                std::cout << "Disease statuses updated" << std::endl;
}

                if (p->control_on == true){
                    Status.updateControl(t); // control and file status updates
if(verbose>1){
                std::cout << "Control statuses updated" << std::endl;
}
                }
                Status.get_premsWithStatus("inf", focalFarms);	// set focalFarms as all farms with disease status inf
                Status.get_premsWithStatus(p->statuses_to_generate_shipments_from, focalFarmsShipments); //Build the set of farms to generate shipments from.

if(verbose>0){
                std::cout <<"Timestep "<<t<<": "
                <<Status.numPremsWithStatus("sus")<<" susceptible, "
                <<Status.numPremsWithStatus("exp")<<" exposed, "
                <<focalFarms.size()<<" infectious, "
                <<Status.numPremsWithStatus("imm")<<" immune premises. "<<std::endl;
                if (p->control_on == true){
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
                Status.newNotSus(notSus); // gets newly not-susceptible farms to remove from consideration
          gridCheck.stepThroughCells(focalFarms,notSus,t); // records new exposures via Status_manager::addPremForEval

                std::clock_t gridcheck_end = std::clock();
                double gridCheckTimeMS = 1000.0 * (gridcheck_end - gridcheck_start) / CLOCKS_PER_SEC;

if(verbose>0){
        		std::cout << "CPU time for checking grid: " << gridCheckTimeMS << "ms." << std::endl;
}


                // determine shipments

                std::vector<Shipment*> fs; // fs = farm shipments, where new shipments are saved
                fs.reserve(10000);
                double shipTimeMS = 0.0;
                if(!focalFarmsShipments.empty() and p->shipments_on) //Only generate shipments if there are any farms to generate from.
                {
                    std::clock_t ship_start = std::clock();
                    std::string time_period = G.get_time_period(); //Current time period we are in given timestep (i.e Q1. Q2, ...)
                    size_t days_in_period = G.get_days_in_period(); //Number of days in this time period.
                    size_t days_rem = G.get_rem_days_of_period(); //Number of days remaining of this time period.
                    Ship.makeShipmentsMultinomial(t, days_in_period, days_rem, time_period, fs,
                                             focalFarmsShipments, G.get_farm_types());
                    std::clock_t ship_end = std::clock();
                    shipTimeMS = 1000.0 * (ship_end - ship_start) / CLOCKS_PER_SEC;
                    // determine which shipments escape ban and are to susceptible farms
                    Status.filter_shipments(fs, t); // Can output shipping from complete "fs" after this point
                    if(verbose>0){std::cout << "CPU time for shipping: " << shipTimeMS << "ms." << std::endl;}
                }

                // evaluate if premises-level control prevents any exposures, expose accordingly
                Status.eval_exposure(t);
                if(verbose>1){std::cout << "SM: eval_exposure complete" << std::endl;}
                Status.add_waitlistMembers(t);
                if(verbose>1){std::cout << "SM: add_waitlistMembers complete" << std::endl;}
				Status.update_ControlResources(t);
                if(verbose>1){std::cout << "SM: update_ControlResources" << std::endl;}
                Status.add_implemented(t);
                if(verbose>1){std::cout << "SM: add_implemented" << std::endl;}

                // at the end of this transmission day, statuses are now...
                Status.get_premsWithStatus("inf", focalFarms); // assign "inf" farms as focalFarms
                int numSuscept = Status.numPremsWithStatus("sus");
                int numExposed = Status.numPremsWithStatus("exp");

                //if bTB like infection is on
                if(p->partial==2){
                    for (auto& f:notSus){ //for the newly exposed farms
                        Pop.verify_currentSize(f); //check if there is a current size, set one if not
                    }
//                    for (auto& f:focalFarms){ //for infected farms
//                        Pop.addAnimals(f, t);
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

                potentialTx = ((focalFarms.size()>0 && numSuscept>0) || (numExposed>0 && numSuscept>0));
                if (p->useMaxPrems==1){
                    int totalInf = Status.get_totalPremsWithStatus("inf");
                    if (totalInf > p->maxInfectiousPrems){
                        potentialTx = 0;
                    }
                }

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

        if (p->printSummary > 0){
            // output summary to file (rep, days inf, run time)
            // rep, # farms infected, # days of infection, seed farm and county, run time
            std::string sumOutFile = batchDateTime;
            sumOutFile += "_summary.txt";
            if (r==1){
                std::string header = "Rep\tNum_Inf\tnAffCounties\tDuration\tSeed_Farms\tSeed_FIPS\tRunTimeSec";
                if (p->control_on == 1){
                for(auto& ct:(p->controlTypes)){
                	std::string ctString = "\t"+ct+"Implemented";
                	ctString += "\t"+ct+"Effective";
                	header+=ctString;
                	for(auto& dct:(p->dcControlTypes)){
                		if (dct==ct){
                			std::string dcImp = "\t"+ct+"ImplementedDCSubset";
											header+=dcImp;
                		}
                	}

                }
                }
                if (p->dangerousContacts_on == 1){
									std::string meandc = "\tmeanDCsPerRP";
									header+=meandc;
								}
                header+="\n";
                printLine(sumOutFile,header);
            }
            std::string repOut = Status.formatRepSummary(r,t,repTimeMS);
            printLine(sumOutFile,repOut);
            }
if(verbose>0){
            std::cout << "Replicate "<< r << " complete." << std::endl<< std::endl;
}
    } // end for loop
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
        G.initShippingParameters(1, p->start_day);
        Shipment_manager Ship(fipsmap, fipsSpeciesMap, p->shipPremAssignment, p->species, p);
        std::vector<std::string> network_out_files;
        std::string netw_number_s = std::to_string(netw_number);
        for(std::string s : p->species)
        {
            network_out_files.push_back(p->batch + "_" + s + "_" + netw_number_s);
        }
        Ship.makeNetwork(network_out_files, G);
    }
}

    std::clock_t process_end = std::clock();
    double process_time = 1000.0 * (process_end - process_start) / CLOCKS_PER_SEC;
    std::cout << "Entire process was alive for " << process_time << "ms." << std::endl;

    return 0;
} // end main()
