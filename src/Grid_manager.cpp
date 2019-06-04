// included in Grid_manager.h: algorithm, queue, stack, unordered_map, vector
#include <cmath> // std::sqrt
#include <fstream>
#include <iostream>
#include <map> // just for printing in order
#include <set> // for comparing otherwise unsorted lists of infected farms, pw vs gridding
#include <sstream>
#include <string>
#include <ctime> // for timing
#include <exception>
#include <algorithm>
// included in Grid_manager.h: grid_cell, farm, shared_functions, tuple, utility
#include "Grid_manager.h"
#include "State.h"
#include "County.h"
#include "shared_functions.h"

/// Loads premises from file, calculates summary statistics
Grid_manager::Grid_manager(const Parameters* p)
	:
	parameters(p),
	speciesOnPrems(p->species),
	susExponents(p->susExponents),
	infExponents(p->infExponents),
	susValues(p->susConsts),
	infValues(p->infConsts),
	kernel(p->kernel),
	committedFarms(0),
	printCellFile(p->printCells),
	batch(p->batch),
	USAMM_temporal_index(-1),
	start_day(p->start_day)
{
	verbose = parameters->verboseLevel;
	FIPS_map.reserve(3200);
	FIPS_vector.reserve(3200);

    // Determine if shipping is turned off in parameters
	shipments_on = p->shipments_on;
	shipment_kernel_str = p->shipment_kernel;
    readFarms(parameters->premFile); // read Farms, make Counties

	if (!shipments_on){
        for(County* c : FIPS_vector){
            // designate as initialized, data not needed if shipping is off
            c -> set_position(0.0, 0.0);
            c -> set_area(0.0);
        }
	}
    readFips_and_states(); // read centroids, areas, state codes, make states

if (verbose>1){
	std::cout << "x min = " << std::get<0>(xylimits) << std::endl;
	std::cout << "x max = " << std::get<1>(xylimits) << std::endl;
	std::cout << "y min = " << std::get<2>(xylimits) << std::endl;
	std::cout << "y max = " << std::get<3>(xylimits) << std::endl;
}

}


Grid_manager::~Grid_manager()
{
  for (auto s : state_map){delete s.second;}
  for (auto c : FIPS_map){delete c.second;}
  for (auto f : farm_map){delete f.second;}
	for (auto gc : allCells){delete gc.second;}
	for (auto ft : farm_types_by_herd){delete ft.second;}
}


void Grid_manager::initShippingParameters(int t, int start_day_in)
{
    //Reads new shipping parameters and updates all state and county specific
    //variables accordingly. Also resets N_todo in states if using usamm v1 parameters.
    start_day = start_day_in;
    usamm_parameters.clear();
    for(size_t i = 0; i < parameters->species.size(); i++)
    {
        std::string this_species = parameters->species.at(i);
        Farm_type* this_ft = farm_types_by_name.at(this_species);
        USAMM_parameters up(parameters, this_ft);
        usamm_parameters[this_ft] = up;
    }
    initFipsCovariatesAndCounties();
    updateShippingParameters(t, true);
}

void Grid_manager::readFips_and_states()
{
    std::clock_t fips_load_start = std::clock();
    int n_counties_loaded = 0;
    int n_states_loaded = 0;

    std::ifstream f(parameters->fipsFile);
    if(f.is_open()) {
	    skipBOM(f);
        while(!f.eof()) {
            std::string line;
			getline(f, line); // get line from file "f", save as "line"
			std::vector<std::string> line_vector = split(line, '\t'); // separate by tabs

            if(! line_vector.empty()) // if line_vector has something in it
            {
                std::string county = line_vector[0];
                std::string state = line_vector[1];
                int int_fips = stringToNum<int>(line_vector[2]);
                if(int_fips == 46113) //Shannon county (46113) changed to Oglala, 46102.
                {
                    int_fips = 46102;
                }
                std::string fips = std::to_string(int_fips);

                int state_code = stringToNum<int>(line_vector[2].substr(0,2));
                double area = stringToNum<double>(line_vector[3]);
                double x = stringToNum<double>(line_vector[4]);
                double y = stringToNum<double>(line_vector[5]);

                if (FIPS_map.find(fips) != FIPS_map.end())  // if fips is present in map
                {
                    // only continue if this county was initiated (because it has 1+ farms)
                    FIPS_map.at(fips) -> set_position(x, y);
                    FIPS_map.at(fips) -> set_area(area);
                    // Add county to its corresponding state object
                    if (state_map.find(state) == state_map.end())  // if state is not present
                    {
                        state_map[state] = new State(state, state_code, parameters->usamm_version);
                        state_vector.push_back(state_map.at(state));
                        n_states_loaded += 1;
                    }
                    if (state_map.find(state) != state_map.end())
                    {
                        FIPS_map.at(fips) -> set_parent_state(state_map.at(state));
                    }
                    else
                    {
                        std::cout << "ERROR: Can't find state that should have just been added: "
                                  << state << std::endl;
                    }
                } // end "if fips is in FIPS_map"
            } // end "if line vector has something in it"
        } // end "while not end of file"
        f.close();
    } else {
        std::cout << "County file not found. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

  n_counties_loaded = FIPS_map.size();

  std::clock_t fips_load_end = std::clock();
  if (verbose > 0){
  std::cout << n_counties_loaded << " counties and " << n_states_loaded << " states created successfully in " <<
              1000.0 * (fips_load_end - fips_load_start) / CLOCKS_PER_SEC <<
              "ms." << std::endl;
  }
}

void Grid_manager::readFarms(const std::string& farm_fname)
{
    //std::clock_t farm_load_start = std::clock();
    farm_map.reserve(850000);
    farm_vector.reserve(850000);
	std::cout<<std::endl;
	int id;
	double x, y;
	std::string fips;
	int fcount = 0;
	std::unordered_map<std::string,double> sumSp; ///> Total count of each species
	std::unordered_map<std::string,double> sumP; ///> Sum of (species count on premises^p), over all premises
	std::unordered_map<std::string,double> sumQ; ///> Sum of (species count on premises^q), over all premises

	// initialize each species sum to 0
	for (std::string s:speciesOnPrems){
		sumSp[s] = 0.0;
		sumP[s] = 0.0;
		sumQ[s] = 0.0;
	}

	std::ifstream f(parameters->premFile);
	if(!f){std::cout << "Premises file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
	if(f.is_open())
	{
	    skipBOM(f);
        if (verbose>0){std::cout << "Premises file open, loading premises." << std::endl;}

		while(! f.eof())
		{
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			std::vector<std::string> line_vector = split(line, '\t'); // separate by tabs

			if(! line_vector.empty()) // if line_vector has something in it
			{
				id = stringToNum<int>(line_vector[0]);
				//Convert string to int, and then back again in order to remove any zeros in the beginning.
				int int_fips = stringToNum<int>(line_vector[1]);
				if(int_fips == 46113) //Shannon county (46113) changed to Oglala, 46102.
                {
                    int_fips = 46102;
                }
				fips = std::to_string(int_fips);

				if (parameters->reverseXY){ // file is formatted as: lat, then long (y, then x)
					y = stringToNum<double>(line_vector[2]);
					x = stringToNum<double>(line_vector[3]);
				} else if (!parameters->reverseXY){ // file is formatted as: long, then lat (x, then y)
					x = stringToNum<double>(line_vector[2]);
					y = stringToNum<double>(line_vector[3]);
				}

				// add species counts - check that number of columns is as expected
				if (line_vector.size() < 4+speciesOnPrems.size()){
					std::cout<<"ERROR (premises file & config 44-46) at line"<<std::endl;
					std::cout<< line <<std::endl;
					std::cout<< ": number of columns with animal populations don't match up with list of species provided."<<std::endl;
					std::cout<<"Exiting..."<<std::endl;
					exit(EXIT_FAILURE);
				}

				unsigned int colcount = 4; // animal populations start at column 4
				unsigned int total_animals = 0;
				std::vector<int> animal_numbers;
				for(size_t i = colcount; i < line_vector.size(); i++)
                {
                    unsigned int tempsize = stringToNum<int>(line_vector[i]);
                    total_animals += tempsize;
                    animal_numbers.push_back(tempsize);
                }

                if(total_animals < 1)
                {
                    if(verbose>1)
                    {
                        std::cout << "Warning: Premises with id " << id << " has no animals of any species. "
                                  << "Skipping this farm..." << std::endl;
                    }
					continue;
                }

				// write farm pointer to private var farm_map
				farm_map[id] = new Farm(id, x, y, fips);
				farm_vector.push_back(farm_map.at(id));
				++fcount;

                std::string herd(speciesOnPrems.size(), '0');
				for (size_t i = 0; i < speciesOnPrems.size(); i++){ // for each species
					std::string sp = speciesOnPrems[i]; //Name of this species
					unsigned int number = animal_numbers[i]; //Number of individuals of this species/type
					farm_map.at(id)->set_speciesCount(sp, number);// set number for species at premises
					sumSp[sp] += number;
					// get infectiousness ("p") for this species
					double p = infExponents.at(sp);
					sumP[sp] += pow(double(number),p);
					// get susceptibility ("q") for this species
					double q = susExponents.at(sp);
					sumQ[sp] += pow(double(number),q);

					// if there are animals of this species, add to fips-species list to sort by population later
					if (number > 0){
						fipsSpeciesMap[fips][sp].emplace_back(farm_map.at(id));
						herd[i] = '1';
					}
				}

				//Assign the correct farm type to the farm.
                Farm_type* farm_type = get_farm_type_by_herd(herd);
                farm_map.at(id)->set_farm_type(farm_type);

// At this point, fixed data pertaining to the farm (coordinates, numbers, types) should
// be final - farm will be copied into other locations (maps by counties)

				// If county doesn't exist, create it
				if (FIPS_map.count(fips) == 0){
                    County* new_county = new County(fips, shipment_kernel_str);
			    	FIPS_map[fips] = new_county;
			    	FIPS_vector.emplace_back(new_county);
                }

				// Add farm to its corresponding county object
				try{
					FIPS_map.at(fips)->County::add_farm(farm_map.at(id));
				}
				catch(std::exception& e){
					std::cout << "When adding premises " << id << " to county " << fips << "." <<
							  e.what() << std::endl;
					exit(EXIT_FAILURE);
 				}


				// compare/replace limits of xy plane
				if (fcount>1){// if this is not the first farm
					if (x < std::get<0>(xylimits)){std::get<0>(xylimits) = x;} // x min
					else if (x > std::get<1>(xylimits)){std::get<1>(xylimits) = x;} // x max

					if (y < std::get<2>(xylimits)){std::get<2>(xylimits) = y;} // y min
					else if (y > std::get<3>(xylimits)){std::get<3>(xylimits) = y;} // y max
					}
				else {
if (verbose>1){std::cout << "Initializing xy limits.";}
					xylimits = std::make_tuple(x,x,y,y);
					// initialize min & max x value, min & max y value
                }

			} // close "if line_vector not empty"
		} // close "while not end of file"
		f.close();
	} // close "if file is open"

    // copy farmlist from farm_map (will be changed as grid is created)
    if (verbose>1){std::cout << "Copying farms from farm_map to farmList..." << std::endl;}
	for (auto& prem: farm_map) {farmList.emplace_back(prem.second);} // "second" value from map is Farm pointer

	// sort farmList by ID for faster matching/subset removal
	std::sort(farmList.begin(),farmList.end(),sortByID<Farm*>);
	// sort within each FIPS map element by farm size (population) for each species
	// (for use in Shipment manager, but this pre-calculates to save run time)
	for (auto& sp:speciesOnPrems){
		std::sort(fipsSpeciesMap[fips][sp].begin(),fipsSpeciesMap[fips][sp].end(),comparePop(sp)); // comparePop struct defined in Grid_manager.h
	}

	// calculate normInf and normSus
	for (auto& sp:speciesOnPrems){
		// sp is species name, sumP is sum of each (herd size^p)
		normInf[sp] = 0.0; //Start at zero and only calc if the total sum of animals of this species is > 0. Otherwise div. by zero and res = nan.
		normSus[sp] = 0.0;
		if(sumSp.at(sp) > 0)
        {
            normInf[sp] = infValues.at(sp)*(sumSp.at(sp)/sumP.at(sp)); // infectiousness normalizer
            normSus[sp] = susValues.at(sp)*(sumSp.at(sp)/sumQ.at(sp)); // susceptibility normalizer
        }


if (verbose>0){
        std::cout<<sp<<" normalized inf: "<<normInf.at(sp)<<", normalized sus: "<<normSus.at(sp)<<std::endl;
}
	}

    double maxFarmInf = 0.0;
    double maxFarmSus = 0.0;
	// calculate and store farm susceptibility and infectiousness
	for (auto& f:farm_map){
		set_FarmSus(f.second);
		set_FarmInf(f.second);
      if (f.second->get_inf_max() > maxFarmInf){maxFarmInf = f.second->get_inf_max();}
    	if (f.second->get_sus_max() > maxFarmSus){maxFarmSus = f.second->get_sus_max();}
	}
if (verbose>0){
    std::cout<<"Max farm inf: "<<maxFarmInf<<", max farm sus: "<<maxFarmSus<<std::endl;
}



}

void Grid_manager::initFipsCovariatesAndCounties()
{
    for(County* c : FIPS_vector)
    {
        c->set_covariates(usamm_parameters);
        c->set_all_counties(&FIPS_vector);
    }
}

void Grid_manager::updateFipsShipping(std::string time_period)
{
    for(County* c : FIPS_vector)
    {
        c->unset_shipping_probabilities();
        c->update_covariate_weights(usamm_parameters, time_period);
    }
    normalizeShippingWeights();
}

void Grid_manager::updateStateLambdas()
{
    for(auto state_pair : state_map)
    {
        for(Farm_type* this_ft : farm_types_vec)
        {
            state_pair.second->update_shipping_rate(this_ft);
        }
    }
}


void Grid_manager::updateStatesShipping(std::string time_period, size_t days_in_period, size_t days_rem,
                                        bool reset)
{
    for(auto state_pair : state_map)
    {
        if(reset)
        {
            state_pair.second->reset_N_todo();
        }
        for(Farm_type* this_ft : farm_types_vec)
        {
            state_pair.second->set_shipping_parameters(parameters->usamm_version, usamm_parameters.at(this_ft),
                                                       this_ft, time_period, days_in_period, days_rem);
        }
    }
}

std::vector<Farm*> Grid_manager::getFarms(std::tuple<int,double,double,double>& cellSpecs, const unsigned int maxFarms/*=0*/)
// based on cell specs, finds farms in cell and saves pointers to farmsInCell
{
if(verbose>1){std::cout << "Getting farms in cell..." << std::endl;
   std::cout << "cellSpecs: " << std::get<0>(cellSpecs) <<", "<< std::get<1>(cellSpecs)
    <<", "<< std::get<2>(cellSpecs) <<", "<< std::get<3>(cellSpecs) << std::endl;}

    // cellSpecs[0] is placeholder for ID number, added when committed
    double x = std::get<1>(cellSpecs);
    double y = std::get<2>(cellSpecs);
    double s = std::get<3>(cellSpecs);
    std::vector<Farm*> inCell;

    // look for farms in cell, those falling on grid boundaries are included, will be removed from list when cell is committed to avoid double counting

    for (auto i:farmList){
    	if ((i->get_x() >= x) && (i->get_x() <= x+s) // if within x bounds of cell
    		&& (i->get_y() >= y) && (i->get_y() <= y+s)) // and within y bounds of cell
    		{ // farm is within the cell
    		inCell.emplace_back(i);
    		if (maxFarms!=0 && inCell.size() > maxFarms){break;}
    		// saves time on retrieving farms if cell will be split & re-checked anyway
    		}
    }
    // (pointers to) farms in inCell should still be sorted by x-coordinate
    return(inCell);
}

void Grid_manager::removeParent(std::stack< std::tuple<int,double,double,double> >& queue)
// The parent cell is the working cell 1st in the queue, so remove first element
{
    queue.pop();
}

void Grid_manager::addOffspring(std::tuple<int,double,double,double> cellSpecs,
	std::stack< std::tuple<int,double,double,double> >& queue)
// offspring cells are quadrants of parent cell
{
	// cellSpecs[0] is placeholder for ID number
    double x = std::get<1>(cellSpecs);
    double y = std::get<2>(cellSpecs);
    double s = std::get<3>(cellSpecs);

    // lower left quadrant: same x/y, side/2
    auto lowerLeft = std::make_tuple(0, x, y, s/2);
    // lower right quadrant: add side/2 to x, same y, side/2
    auto lowerRight = std::make_tuple(0, x+s/2, y, s/2);
    // upper left quadrant: same x, add side/2 to y, side/2
    auto upperLeft = std::make_tuple(0, x, y+s/2, s/2);
    // upper right quadrant: add side/2 to x and y, side/2
    auto upperRight = std::make_tuple(0, x+s/2, y+s/2, s/2);

    // add offspring cells to queue (in reverse order so lower left is first)
    queue.emplace(upperRight);
    queue.emplace(upperLeft);
    queue.emplace(lowerRight);
    queue.emplace(lowerLeft);
}

void Grid_manager::commitCell(std::tuple<int,double,double,double> cellSpecs, std::vector<Farm*>& farmsInCell)
// write cellSpecs as class Grid_cell into set allCells
{
	int id;
    double x, y, s;
    id = std::get<0>(cellSpecs);
    x = std::get<1>(cellSpecs);
    y = std::get<2>(cellSpecs);
    s = std::get<3>(cellSpecs);
    std::vector<Farm*> farms = farmsInCell;

    Grid_cell* cellToAdd = new Grid_cell(id, x, y, s, farms);

    allCells.emplace(id,cellToAdd); // add to map of all committed cells, with id as key

    std::set<std::string> countiesInCell = cellToAdd -> Grid_cell::get_counties();
    // store cells associated with each county
    for (auto& c: countiesInCell){
			cellsByCounty[c].emplace_back(cellToAdd);
		}

    committedFarms += farmsInCell.size();
    assignCellIDtoFarms(id,farmsInCell);
    removeFarmSubset(farmsInCell, farmList);
if (verbose>1){std::cout<<"Cell committed, id "<<id<<std::endl;}
}

void Grid_manager::splitCell(std::tuple<int,double,double,double>& cellSpecs,
	std::stack< std::tuple<int,double,double,double> >& queue)
{
    removeParent(queue);
    addOffspring(cellSpecs,queue);
}

void Grid_manager::assignCellIDtoFarms(int cellID, std::vector<Farm*>& farmsInCell)
{
	for (auto& f:farmsInCell){
		f->set_cellID(cellID);
	}
}

void Grid_manager::initiateGrid(const unsigned int in_maxFarms, const int minCutoff)
// maxFarms: If cell contains at least this many farms, subdivision will continue
// minCutoff: minimum cell size
{
	set_maxFarms(in_maxFarms);
	if (verbose > 0){
	std::cout << "Max farms set to " << maxFarms << std::endl;
	}
    if(verbose>0){std::cout << "Splitting into grid cells..." << std::endl;}
 	int cellCount = 0;
    std::stack<std::tuple<int,double,double,double>> queue;// temporary list of cells to check for meeting criteria for commitment
    std::vector<Farm*> farmsInCell; // vector of (pointers to) farms in working cell - using vector to get to specific elements

    double min_x = std::get<0>(xylimits)-0.1;
    double max_x = std::get<1>(xylimits)+0.1;
    double min_y = std::get<2>(xylimits)-0.1;
    double max_y = std::get<3>(xylimits)+0.1;

    double side_x = max_x - min_x;
    double side_y = max_y - min_y;
    if(verbose>1)
    {
    	std::cout << "side_x = " << side_x << std::endl;
    	std::cout << "side_y = " << side_y << std::endl;
    }

    // use whichever diff is larger, x or y
    if (side_y > side_x)
       side_x = side_y;
    if(verbose>1){std::cout << "Using larger value " << side_x << std::endl;}

    // add cell specifications to temporary tuple
    std::tuple<int,double,double,double> cellSpecs = std::make_tuple(cellCount, min_x, min_y, side_x);
    if(verbose>1){std::cout << "cellSpecs: " << std::get<0>(cellSpecs) <<", "<< std::get<1>(cellSpecs)
    	<<", "<< std::get<2>(cellSpecs) <<", "<< std::get<3>(cellSpecs) << std::endl;}

    // add initial cell to the queue
    queue.emplace(cellSpecs);

    // while there are any items in queue
    while(queue.size()>0)
    {
    if(verbose>1){std::cout << std::endl << "Queue length = " << queue.size() << std::endl;}
    cellSpecs = queue.top(); // set first in queue as working cell

	if(verbose>1){
    	std::cout << "Cell side length = " << std::get<3>(cellSpecs) << ". ";
    }

	// Case A: side length of cell is smaller than kernel - immediate commit
	if (std::get<3>(cellSpecs) < minCutoff){ // if side < kernel diameter
		farmsInCell = getFarms(cellSpecs); // want ALL farms, so don't include maxFarms as argument
		if (farmsInCell.size() > 0){ // if there are farms in cell, commit
        	std::get<0>(cellSpecs) = cellCount;
        	commitCell(cellSpecs,farmsInCell);
        	cellCount = cellCount+1;
        	if (verbose==2){std::cout << "Side smaller than kernel diameter.";}
        	queue.pop(); // remove parent cell from front of queue
        } else { // no farms in cell, remove from queue w/o committing
        	queue.pop(); // remove parent cell from front of queue
            if(verbose==2){std::cout << "No farms, removed cell, queue length = " << queue.size() << std::endl;}
        }
    // Case B: side length of cell >= minimum, check farm density and split if needed
    } else if (std::get<3>(cellSpecs) >= minCutoff){ // side >= kernel diameter
    	if(verbose==2){std::cout << "Side bigger than kernel, stepping in..." << std::endl;}
    	farmsInCell = getFarms(cellSpecs, maxFarms); // copy up to maxFarms farms in cell to farmsInCell)
    	if(verbose==2){std::cout << "Farms in cell = " << farmsInCell.size() << std::endl;}
        if (farmsInCell.size() >= maxFarms){
        // if farm density too high, split
			if(verbose==2){std::cout << "Too many farms, splitting cell..." << std::endl;}
			splitCell(cellSpecs,queue);
        }
        else if (farmsInCell.size() > 0 && farmsInCell.size() < maxFarms){
        // farm density is below maximum, commit
            	std::get<0>(cellSpecs) = cellCount;
                commitCell(cellSpecs,farmsInCell);
                cellCount = cellCount+1;
                if (verbose>1){std::cout << "Cell committed: #" << cellCount<<std::endl;}
                queue.pop(); // remove parent cell from front of queue
            }
        else if (farmsInCell.empty()){
        // cell has no farms at all - remove from queue w/o committing
            queue.pop(); // remove parent cell from front of queue
            if(verbose==2){std::cout << "No farms, removed cell, queue length = "
            	<< queue.size() << std::endl;}
        }
    }
    } // end "while anything in queue"

if (verbose){
	std::cout << "Grid of "<< allCells.size()<<" cells created, with min side "<<minCutoff<<
	" and max "<<maxFarms<<" farms. Pre-calculating distances..." << std::endl;
}
	if (farm_map.size()!=committedFarms){
		std::cout<<"ERROR: "<<committedFarms<<" farms committed, expected "
		<<farm_map.size()<<std::endl;
		for (auto& dropped:farmList){std::cout<<"ID "<<dropped->get_id()<<", x "
			<<dropped->get_x()<<", y "<<dropped->get_y()<<std::endl;}
		exit(EXIT_FAILURE);
		}
	makeCellRefs();
	if(verbose>0){std::cout << "Grid initiated using density parameters. ";}
	if (printCellFile > 0){printCells();}
}

void Grid_manager::initiateGrid(std::string& cname)
// overloaded (alternate) constructor that reads in external file of cells
{
	// read in file of premises
	std::vector<Farm*> farmsInCell;

	std::ifstream f(cname);
	if(!f){std::cout << "Input file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
	if(f.is_open())
	{
	if (verbose>0){std::cout << "File open" << std::endl;}
		while(! f.eof())
		{
			std::tuple<int,double,double,double> cellSpecs;
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			std::vector<std::string> line_vector = split(line, '\t'); // separated by tabs

			if(! line_vector.empty()) // if line_vector has something in it
			{ // convert each string piece to double
				if (verbose==2){std::cout << "Reading cell: ";}
				std::get<0>(cellSpecs)=stringToNum<int>(line_vector[0]); //id
				if (verbose==2){std::cout << std::get<0>(cellSpecs) << ", ";}
				std::get<1>(cellSpecs)=stringToNum<double>(line_vector[1]); //x
				if (verbose==2){std::cout << std::get<1>(cellSpecs) << ", ";}
				std::get<2>(cellSpecs)=stringToNum<double>(line_vector[2]); //y
				if (verbose==2){std::cout << std::get<2>(cellSpecs) << ", ";}
				std::get<3>(cellSpecs)=stringToNum<double>(line_vector[3]); //side
				if (verbose==2){std::cout << std::get<3>(cellSpecs) << ". ";}
				// line_vector[4] is num farms-ignored (gets reassigned)
				farmsInCell = getFarms(cellSpecs);
				if (verbose==2){std::cout << farmsInCell.size() << " farms assigned to cell." <<
					std::endl;}
				if(farmsInCell.empty()){
					std::cout << "Cell " << std::get<0>(cellSpecs) << " has no farms - ignoring." <<
						std::endl;
					// cell will not be added to list
					}
				else if (!farmsInCell.empty()){
				// save cell with farms within
					allCells[std::get<0>(cellSpecs)] = new Grid_cell(std::get<0>(cellSpecs),
						std::get<1>(cellSpecs), std::get<2>(cellSpecs), std::get<3>(cellSpecs), farmsInCell);
					assignCellIDtoFarms(std::get<0>(cellSpecs),farmsInCell);
					removeFarmSubset(farmsInCell, farmList);
					}
			} // close "if line_vector not empty"
		} // close "while not end of file"
	} // close "if file is open"
	f.close();
	if (verbose>0){std::cout << "File closed" << std::endl;}
	std::cout << allCells.size() << " cells loaded from file."<<std::endl;
	if (!farmList.empty()){
		Farm* f = farmList[0];
		std::cout << farmList.size() << " unassigned farms, first: " << f->get_id() << ": x=" << f->get_x() <<
			", y=" << f->get_y() << std::endl;
		}
	makeCellRefs();
	if (printCellFile > 0){printCells();}
}

void Grid_manager::initiateGrid(double cellSide)
{
    double min_x = std::get<0>(xylimits);
    double max_x = std::get<1>(xylimits);
    double min_y = std::get<2>(xylimits);
    double max_y = std::get<3>(xylimits);

    std::unordered_map<int, std::vector<Farm*>> cellFarmMap;
    std::vector<double> xlist, ylist; // list of each x corner, y corner
	std::vector<int> uniquex; // list of elements of first unique x values
		uniquex.emplace_back(0); // include first value (element/index 0)
	int cellCount = 0;
    // all x points will be from min_x to max_x by cellSide
    // let max be max + cellside for extra wiggle room, cells w/o farms will be excluded later
    for (auto x = min_x; x <= (max_x+cellSide); x+=cellSide)
    {
    	// all y points will be from min_y to max_y by cellSide
    	for (auto y = min_y; y <= (max_y+cellSide); y+=cellSide)
    	{
    		xlist.emplace_back(x);
    		ylist.emplace_back(y);
    		++cellCount;
        	}
    uniquex.emplace_back(cellCount); // new x-value will start at element cellCount
    }
    // remove last x value (was added on after all finished)
    uniquex.pop_back();
    // assign farms to cells (have whole list already, so don't use getFarms)
    // compare each farm's coordinates to:
    // unique cell x values (increment according to xChanges)
    // y values once x value is found

    std::vector<Farm*> farmListByX(farmList); // farmListByX is a copy of farmList
    std::sort(farmListByX.begin(), farmListByX.end(), sortByX);
    // indices for moving around the cell list
    int xi = 0;
    int i = 0;
    int fcount = 0;

    std::vector<int> seedFarms;

    for (auto& f:farmListByX)
    {
		double farmx = f->Farm::get_x();
		double farmy = f->Farm::get_y();
    	bool cellFound = 0;
    	fcount++;

    	while(!cellFound){
    		// if farm x is...
			if (farmx >= xlist[uniquex[xi]] && farmx < xlist[uniquex[xi+1]]) // between this and next unique x value
			{ // move ahead and narrow down y value
				i = uniquex[xi];
				while (!cellFound){
				if (farmy >= ylist[i] && farmy < ylist[i+1])
				// within y range - this is the cell
				{
					cellFarmMap[i].emplace_back(f);
					cellFound = 1;
				}
				else if (farmy >= ylist[i]){++i;}
				else if (farmy < ylist[i]){std::cout << "Farm in previous x value range, "; --i;}
				else {std::cout << "Farm y range not found." << std::endl;}
				} // end 2nd while cell not found
			}
			else if (farmx >= xlist[uniquex[xi+1]]) // or if farm x >= next larger x value, increase xi
			{++xi;}
			else if (farmx < xlist[i]) // or if farm x < this x value
			{std::cout << "Farm in previous x value range, "; --xi;}
			else {std::cout << "Farm x range not found." << std::endl;}
		} // end 1st while cell not found
    } // end for each farm
    std::cout << "Done placing farms." << std::endl;

   // commit all cells with farms
   bool printNumFarms = 0;
   std::string allLinesToPrint;
   int actualCellCount = 0;
   for (auto c=0; c!=cellCount; ++c)
   {
	   if (cellFarmMap[c].size()>0){ // if there are any farms in this cell:
	   		allCells[actualCellCount] = new Grid_cell(actualCellCount, xlist[c], ylist[c], cellSide, cellFarmMap[c]);
	   		assignCellIDtoFarms(actualCellCount,cellFarmMap[c]);
	   		++actualCellCount;
	   }
   }
   if (printNumFarms){
   		std::string ofilename = "farmsPerUnifCell.txt";
   		std::ofstream f(ofilename);
	if(f.is_open())
	{
		f << allLinesToPrint;
		f.close();
	}
	}

	std::cout << "Grid loaded with " << actualCellCount << " uniform cells. Pre-calculating distances..." << std::endl;
	makeCellRefs();
	if (printCellFile > 0){printCells();}
}

std::string Grid_manager::to_string(Grid_cell& gc) const
// overloaded to_string function, makes tab-delim string (one line) specifically for cell
{
	std::string toPrint;
	char temp[20];
	std::vector<double> vars;
		vars.resize(5);
		vars[0] = gc.get_id();
		vars[1] = gc.get_x();
		vars[2] = gc.get_y();
		vars[3] = gc.get_s();
		vars[4] = gc.get_num_farms();

	for(auto it = vars.begin(); it != vars.end(); ++it)
	{
		sprintf(temp, "%f\t", *it);
		toPrint += temp;
	}

	toPrint.replace(toPrint.end()-1, toPrint.end(), "\n");

	return toPrint;
}

/// Used in commitCell in grid initiation.
void Grid_manager::removeFarmSubset(std::vector<Farm*>& subVec, std::vector<Farm*>& fullVec)
{
	unsigned int expectedSize = fullVec.size()-subVec.size();
//	std::cout << "Removing "<<subVec.size()<<" farms from list of "<<fullVec.size()<<std::endl;

	// put vectors into fips-indexed maps to speed up matching
	std::unordered_map< std::string, std::vector<Farm*> > subMap, fullMap;
	for (auto& sv:subVec){
		subMap[sv->get_fips()].emplace_back(sv);}
	for (auto& fv:fullVec){
		fullMap[fv->get_fips()].emplace_back(fv);}

	for (auto& sub:subMap){
		// for each fips in subset list
		std::string fips = sub.first;
		// if needed, sort both lists of farms in this FIPS, by ID
		std::sort(sub.second.begin(),sub.second.end(),sortByID<Farm*>);
		std::sort(fullMap.at(fips).begin(),fullMap.at(fips).end(),sortByID<Farm*>);
		// iterate through full list, erasing matching sub as found
		auto it2 = fullMap.at(fips).begin();
		for(auto it = sub.second.begin(); it != sub.second.end(); it++){
		// loop through each farm in this FIPS
			while (it2 != fullMap.at(fips).end()){ // while end of full list not reached
				if(*it2 == *it){ // finds match in farmList to farmInCell
					fullMap.at(fips).erase(it2); // remove from farmList
					break; // start at next farm instead of looping over again
				}
				it2++;
			}
		}
	}
	// rewrite fullVec
	std::vector<Farm*> temp;
	for (auto& f1:fullMap){
	  for (auto& f2:f1.second){
		temp.emplace_back(f2);}}
	fullVec = temp;

	if (expectedSize != fullVec.size()){
		std::cout << "Error in removeFarmSubset: expected size"<< expectedSize <<
		", actual size: "<< fullVec.size() <<". Exiting...";
		exit(EXIT_FAILURE);
	}

}

/// Gets shortest distance between two cells and squares that distance (for input to kernel)
double Grid_manager::shortestCellDist2(Grid_cell* cell1, Grid_cell* cell2)
{
	double cellDist2 = 0; // squared distance between cells
 	int cell1_id = cell1->Grid_cell::get_id();
 	int cell2_id = cell2->Grid_cell::get_id();
  if (cell1_id != cell2_id){ // else if comparing to self, distance is already 0
	double cell1_x = 0;
	double cell1_y = 0;
	double cell2_x = 0;
	double cell2_y = 0; // will use these points to calc distance

	double cell1_South = cell1->Grid_cell::get_y(); // lower boundary of cell1
	double cell1_North = cell1->Grid_cell::get_y()+cell1->Grid_cell::get_s(); // upper boundary of cell1
	double cell1_West = cell1->Grid_cell::get_x(); // leftmost boundary of cell1
	double cell1_East = cell1->Grid_cell::get_x()+cell1->Grid_cell::get_s(); // rightmost boundary of cell1

	double cell2_South = cell2->Grid_cell::get_y(); // lower boundary of cell2
	double cell2_North = cell2->Grid_cell::get_y()+cell2->Grid_cell::get_s(); // upper boundary of cell2
	double cell2_West = cell2->Grid_cell::get_x(); // leftmost boundary of cell2
	double cell2_East = cell2->Grid_cell::get_x()+cell2->Grid_cell::get_s(); // rightmost boundary of cell2

	// In comparing cell positions, due to nestedness of sub-cells, cell 1 could be:
	// Horizontally: W of, E of, or directly above/below all or part of cell2.
	// Vertically: N of, S of, or directly beside all or part of cell2.

	// Determine horizontal relationship and set x values accordingly:
	if(verbose>1){std::cout << "Cell " << cell1_id << " is ";}

	if (cell1_East <= cell2_West) // cell1 west of cell2
		{
 		if(verbose>1){std::cout << "west of and ";}
		cell1_x = cell1_East;
		cell2_x = cell2_West;
		}
	// or cell1 is east of cell2
	else if (cell1_West >= cell2_East)
		{
 		if(verbose>1){std::cout << "east of and ";}
		cell1_x = cell1_West;
		cell2_x = cell1_East;
		}
	// or cell1 is directly atop all or part of cell2
	else // if ((cell1_East > cell2_West) && (cell1_West < cell2_East))
		{
 		if(verbose>1){std::cout << "vertically aligned with and ";}
		//cell1_x = 0; // already initialized as 0
		//cell2_x = 0; // already initialized as 0
		// only use distance between y values
		}

	// Determine vertical relationship and set y values accordingly:
	if (cell1_South >= cell2_North) // cell1 north of cell2
		{
 		if(verbose>1){std::cout << "north of cell "<< cell2_id << std::endl;}
		cell1_y = cell1_South;
		cell2_y = cell2_North;
		}
	// or cell1 is below cell2
	else if (cell1_North <= cell2_South)
		{
 		if(verbose>1){std::cout << "south of cell "<< cell2_id << std::endl;}
		cell1_y = cell1_North;
		cell2_y = cell2_South;
		}
	// or cell1 is directly beside cell2
 	else // if ((cell1_South < cell2_North) && (cell1_North > cell2_South))
		{
 		if(verbose>1){std::cout << "horizontally aligned with cell "<< cell2_id << std::endl;}
		//cell1_y = 0; // already initialized as 0
		//cell2_y = 0; // already initialized as 0
		// only use distance between x values
		}

	double xDiff = cell1_x-cell2_x;
	double yDiff = cell1_y-cell2_y;
	cellDist2 = xDiff*xDiff + yDiff*yDiff;

  } // end if cells 1 and 2 are different

return cellDist2;
}

/// Calculates kernel values * max susceptibility (susxKern) for each pair of Grid_cells.
/// Stored with each Grid_cell is a map with all other cells as keys, with values susxKern
/// Neighbors (other Grid_cells with shortest distance = 0) are also stored with each Grid_cell
void Grid_manager::makeCellRefs()
// Although all the ID referencing seems a bit much, this is one way to ensure the order of the cells checked
{
	std::unordered_map<Grid_cell*, std::unordered_map<int, double> > susxKern;

	for (unsigned int whichCell1=0; whichCell1 != allCells.size(); ++whichCell1){
		Grid_cell* cell1 = allCells.at(whichCell1);
		for (unsigned int whichCell2 = whichCell1; whichCell2 != allCells.size(); ++whichCell2){
			Grid_cell* cell2 = allCells.at(whichCell2);
			// get distance between grid cells 1 and 2...
			// if comparing to self, distance=0
			double shortestDist2 = 0;
			if (whichCell2 != whichCell1) { // overwrite if cells are different
				shortestDist2 = shortestCellDist2(cell1, cell2);
if(verbose>1){std::cout << "Distance squared between "<<whichCell1<<" & "<<whichCell2<<": "<<shortestDist2<<std::endl;}
			}
			// save adjacent neighbors, not including self (for distance-based control)
			if (shortestDist2 == 0 && whichCell1 != whichCell2){
				cell1->addNeighbor(cell2);
				cell2->addNeighbor(cell1);
			}
			// kernel value between c1, c2
			double gridValue = kernel->atDistSq(shortestDist2);
if(verbose>1){std::cout << "Kernel between "<<whichCell1<<"&"<<whichCell2<<": "<<gridValue<<std::endl;}
				// store kernel * max sus (part of all prob calculations)
				double maxS2 = cell2->Grid_cell::get_maxSus();
				if (parameters->dangerousContacts_on){
					susxKern[cell1][whichCell2] = maxS2 * gridValue * (parameters->maxDCScale);
				} else {
					susxKern[cell1][whichCell2] = maxS2 * gridValue;
				}
if(verbose>1){std::cout << "Stored in-range sus*kernel: cells "<<whichCell1<<" & "<<whichCell2<<", susxKern: "<<
	susxKern.at(cell1).at(whichCell2)<<std::endl;}

				// if not comparing to self, calc/store other direction (this was a big bug - double counting self as neighbor)
				if (whichCell1 != whichCell2){
					double maxS1 = cell1->Grid_cell::get_maxSus();
					if (parameters->dangerousContacts_on){
						susxKern[cell2][whichCell1] = maxS1 * gridValue * (parameters->maxDCScale);
					} else {
 						susxKern[cell2][whichCell1] = maxS1 * gridValue;
 					}
 if(verbose>1){std::cout << "Stored in-range sus*kernel: cells "<<whichCell2<<" & "<<whichCell1<<", susxKern: "<<
 	susxKern.at(cell2).at(whichCell1)<<std::endl;}
				}
		} // end for each cell2
	} // end for each cell1

	// assign kernel maps to individual cells
	for (auto& k:susxKern){
		k.first->take_KernelValues(k.second);
	}


if (verbose>0){std::cout << "Kernel distances and neighbors recorded." << std::endl;}
}

/// Used after grid creation to assign susceptibility values to individual premises
void Grid_manager::set_FarmSus(Farm* f)
{
	// calculates species-specific susceptibility for a premises
	// USDOSv1 uses scaling factor (for q, susceptibility)
	// 2.086 x 10^-7, or that times sum of all US cattle: 19.619
	double premSus = 0.0;
	for (auto& sp:speciesOnPrems){
		double n_animals = double(f->get_size(sp)); // i.e. get_size("beef") gets # of beef cattle on premises
		double spSus = normSus.at(sp)*pow(n_animals, susExponents.at(sp)); // multiply by stored susceptibility value for this species/type
		premSus += spSus; // add this species to the total for this premises
	}
	f->set_sus(premSus);
	//std::cout << "Prem sus = " << f->get_sus() << ". ";
}

/// Used after grid creation to assign infectiousness values to individual premises
void Grid_manager::set_FarmInf(Farm* f)
{
	// calculates species-specific infectiousness for a premises
	// USDOSv1 uses scaling factor (for p, transmissibility)
	// 2.177 x 10^-7, or that times sum of all US cattle: 20.483
	double premInf = 0.0;
	for (auto& sp:speciesOnPrems){
		double n_animals = double(f->get_size(sp)); // i.e. get_size("beef") gets # of beef cattle on premises
		double spInf = normInf.at(sp)*pow(n_animals, infExponents.at(sp)); // susceptibility value for this species/type
		premInf += spInf; // add this species to the total for this premises
	}
	f->set_inf(premInf);
	//std::cout << "Prem inf = " << f->get_inf() << ". ";

}

/// Prints file with specifications of cells
void Grid_manager::printCells()
{
	std::string sumOutFile = batch;
	sumOutFile += "_cells.txt";

	std::string header = "Cell\tLow_X\tLow_Y\tSide_meters\tNum_Farms\tFIPS\n";
	printLine(sumOutFile,header);

	for (auto& c:allCells){
		std::string cellOut;
		addItemTab(cellOut, c.second->Grid_cell::get_id());
		addItemTab(cellOut, c.second->Grid_cell::get_x());
		addItemTab(cellOut, c.second->Grid_cell::get_y());
		addItemTab(cellOut, c.second->Grid_cell::get_s());
		addItemTab(cellOut, c.second->Grid_cell::get_num_farms());

		std::vector<Farm*> inCell = c.second->Grid_cell::get_farms();
		std::set<std::string> counties;
		for (auto& ic:inCell){
			counties.emplace(ic->Farm::get_fips());
		}

		std::string countycomma;
		for (auto& co:counties){
			countycomma += co;
			countycomma += ",";
		}
		countycomma.pop_back();  // remove last comma

		addItemTab(cellOut,countycomma);
		cellOut.replace(cellOut.end()-1, cellOut.end(), "\n"); // add line break at end

		printLine(sumOutFile,cellOut);
	}

}

void Grid_manager::updateShippingParameters(int t, bool restart)
{
    //Get the day of the year given the current timestep and the start day (from config file).
    //This is needed since there can be (probably is) different shipping parameter sets for different
    //time periods (eg. quarters, months) and the transitions need to be kept track of.

    int day_of_year = get_day_of_year(t, start_day);
//    int day_of_year = ((t-2 + start_day) % 365) + 1;
    int current_year = ((t-2 + start_day) / 365) + 1;
    bool new_year = false;
    if(current_year > USAMM_current_year)
    {
        USAMM_current_year = current_year;
        new_year = true;
    }
    //This vector contains the start day for the time periods.
    const std::vector<int>& start_points = parameters->USAMM_temporal_start_points;
    bool update_parameters = false;
    //Update the index (USAMM_temporal_index) that tells us what time period we are in
    //given the current timestep - check each start point to see which period  the current day of year falls into.
    for(size_t i = start_points.size()-1; i >= 0; i--)
    {
        if(day_of_year >= start_points.at(i))
        {
            //If the index given by timestep has changed compared to USAMM_temporal_index we are in a
            //new time period.
            if(i != USAMM_temporal_index or new_year)
            {
                //Time period has changed from previous day. Parameters need to be updated. This will
                //always happen upon construction of Grid_manager since USAMM_temporal_index = -1 on construction.
                //This causes all states to update the relevant parameters and then the counties to
                //update their covariate weights (which are functions of time-dependent covariate parameters)
                //and shipping probabilities.
                USAMM_temporal_index = i;
                USAMM_temporal_name = parameters->USAMM_temporal_order.at(USAMM_temporal_index);
                days_in_period = parameters->USAMM_temporal_n_days.at(USAMM_temporal_index);
                update_parameters = true;
            }
            break;
        }
    }

    //Calculate the number of remaining days in this time period. This is used when determining
    //the number of shipments from a state when starting the simulation/generation in the middle
    //of a time period.
    if(USAMM_temporal_index == start_points.size()-1)
    {
        days_rem_of_period = 366 - day_of_year;
    }
    else
    {
        days_rem_of_period = start_points[USAMM_temporal_index+1] - day_of_year;
    }
    days_rem_of_period = std::min(days_rem_of_period, size_t(parameters->timesteps + 1 - t));

    if(update_parameters or restart)
    {
        updateStatesShipping(USAMM_temporal_name, days_in_period, days_rem_of_period, restart);
        updateFipsShipping(USAMM_temporal_name);
        if(parameters->usamm_version == 2)
        {
            updateStateLambdas();
        }
        else if(parameters->usamm_version < 1 or parameters->usamm_version > 2)
        {
            std::cout << "Shipments turned off or USAMM version undefined, but still in "
                      << "shipment-related code (updateShippingParameters). This should not happen."
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

std::string Grid_manager::get_generation_string(Farm_type* ft)
{
    return usamm_parameters.at(ft).get_generation_string();
}

void Grid_manager::normalizeShippingWeights()
{
    for(Farm_type* current_ft : farm_types_vec)
    {
        for(State* s : state_vector)
        {
            s->normalize_shipping_weights(current_ft);
        }
    }
}

Farm_type* Grid_manager::get_farm_type_by_herd(std::string herd)
{
    if(farm_types_by_herd.find(herd) == farm_types_by_herd.end())
    {
        //This farm type has not been created yet
        unsigned int farm_type_index = farm_types_by_herd.size();
        farm_types_by_herd[herd] = new Farm_type(herd, parameters->species, farm_type_index);
        farm_types_vec.push_back(farm_types_by_herd[herd]);
        farm_types_by_name[farm_types_by_herd.at(herd)->get_species()] = farm_types_by_herd.at(herd);
        return farm_types_by_herd.at(herd);
    }
    else
    {
        return farm_types_by_herd.at(herd);
    }
}

Farm_type* Grid_manager::get_farm_type_by_name(std::string name)
{
    return farm_types_by_name.at(name);
}

std::vector<Farm_type*> Grid_manager::get_farm_types()
{
    return farm_types_vec;
}

std::string Grid_manager::get_time_period()
{
    return USAMM_temporal_name;
}

size_t Grid_manager::get_days_in_period()
{
    return days_in_period;
}

size_t Grid_manager::get_rem_days_of_period()
{
    return days_rem_of_period;
}

/// Reads seed infection source file, expects vector of strings
void Grid_manager::read_seedSource(std::string seedSource, std::vector<std::string>& output)
{
	std::vector<std::string> tempOutput;
	std::ifstream f(seedSource);
	if(!f){std::cout << "Seed source file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
	if(verbose>1){std::cout << "Loading seed data from "<<seedSource<<std::endl;}
		while(! f.eof()){
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			if(! line.empty()){ // if line has something in it
				// check for 5 characters? no straightforward "nchar" function
				tempOutput.emplace_back(line);
			} // close "if line_vector not empty"
		} // close "while not end of file"
		tempOutput.swap(output);
		if(verbose>0){std::cout << " Closed seed file." << std::endl;}
}

/// Overloaded version of readSource, expects vector of ints back
void Grid_manager::read_seedSource(std::string seedSource, std::vector<int>& output)
{

	std::vector<int> tempOutput;
	std::ifstream f(seedSource);
	if(!f){std::cout << "Seed source file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
if(verbose>0){std::cout << "Loading seed data from "<<seedSource<<std::endl;}
		while(! f.eof()){
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			if(! line.empty()){ // if line has something in it
				tempOutput.emplace_back(stringToNum<int>(line));
			} // close "if line_vector not empty"
		} // close "while not end of file"
		tempOutput.swap(output);
if(verbose>0){std::cout << " Closed seed file." << std::endl;}
}

/// Overloaded version of readSource, expects vector of vector of ints back
void Grid_manager::read_seedSource(std::string seedSource, std::vector<std::vector<int>>& output)
{
	std::vector<int> tempLine;
	std::vector<std::vector<int>> tempOutput;
	std::ifstream f(seedSource);
	if(!f){std::cout << "Seed source file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
	if(verbose>0){std::cout << "Loading seed data from "<<seedSource<<std::endl;}
		while(! f.eof()){
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			if(! line.empty()){ // if line has something in it
				tempLine = stringToIntVec(line);
				tempOutput.emplace_back(tempLine);
			} // close "if line_vector not empty"
		} // close "while not end of file"
		tempOutput.swap(output);
		if(verbose>0){std::cout << " Closed seed file." << std::endl;}
}

///< Selects one random premises per county (returns one for every county with at least one
/// premises)
void Grid_manager::select_randomPremisesPerCounty(std::vector<std::vector<Farm*>>& output)
{
if (verbose>0){
	std::cout << "Selecting a random farm from each of "<<FIPS_vector.size()<<" counties."
	<< std::endl;
}
	std::vector<std::vector<Farm*>> tempOutput;
		tempOutput.reserve(FIPS_vector.size());
	for (auto& c:FIPS_vector){
		std::vector<Farm*> premisesInCounty = c->get_farms();
		std::vector<Farm*> tempPremVector;
		random_unique(premisesInCounty, 1, tempPremVector); // stores random premises in tempPremVector
		tempOutput.emplace_back(tempPremVector);
	}
//std::cout<<"Farm "<<tempOutput.front().front()->Farm::get_id()<<" selected for exposure."<<std::endl;
	tempOutput.swap(output);
}

///< Overloaded version of select random premises per county, for specific list of FIPS
void Grid_manager::select_randomPremisesPerCounty(std::vector<std::string> fipsStrings,
	std::vector<std::vector<Farm*>>& output)
{
if (verbose>0){
	std::cout << "Selecting a random farm from each of "<<fipsStrings.size()<<" counties."
	<< std::endl;
}
	std::vector<std::vector<Farm*>> tempOutput;
		tempOutput.reserve(fipsStrings.size());
	std::vector<Farm*> tempPremVector(1);
	for (auto& c:fipsStrings){
        //Convert string to int, and then back again in order to remove any zeros in the beginning.
        int int_fips = stringToNum<int>(c);
        std::string fips = std::to_string(int_fips);
		std::vector<Farm*> premisesInCounty = FIPS_map.at(fips)->get_farms();
		random_unique(premisesInCounty, 1, tempPremVector); // stores random premises in tempPremVector
		tempOutput.emplace_back(tempPremVector);
	}
	tempOutput.swap(output);
}

/// Each run ( = simulation = replicate) is initiated by designating one of the following
/// as "exposed":
/// - one randomly chosen farm per county (seedSource = "allFips" or seedSourceType = "fips")
/// - one specific premises (seedSource = "singlePremises", one ID per line of file)
/// - multiple specific premises (seedSource = "multiplePremises", multiple comma-separated
///   premises IDs per line of file)
void Grid_manager::get_seedPremises(std::vector<std::vector<Farm*>>& output){
	std::vector<std::vector<Farm*>> seedFarmsByRun;
	if (parameters->seedSource == "allFips"){// if seeding one per county
		// get a random farm per county
        select_randomPremisesPerCounty(seedFarmsByRun); // saves to seedFarmsByRun
    } else if (parameters->seedSource != "allFips"){ // will need to read file
    	if (parameters->seedSourceType == "fips"){
    		std::vector<std::string> sourceFipsStrings;
    		// read source file (expect vector of strings back)
    		read_seedSource(parameters->seedSource, sourceFipsStrings);
    		// get a random farm per county
    		select_randomPremisesPerCounty(sourceFipsStrings, seedFarmsByRun); // saves to seedFarmsByRun
    	} else if (parameters->seedSourceType == "singlePremises"){
    		std::vector<int> sourcePremIDs;
    		// read source file (expect vector of ints back)
    		read_seedSource(parameters->seedSource, sourcePremIDs);
			// translate strings to premIDs
			for (auto& pID:sourcePremIDs){
				if (farm_map.count(pID) == 0){
					std::cout << "Warning: Premises " << pID <<
					" not found, skipping this source of infection." << std::endl;
				} else {
					std::vector<Farm*> tempPremVector;
					tempPremVector.emplace_back(farm_map.at(pID));
					seedFarmsByRun.emplace_back(tempPremVector);
				}
			}
    	} else if (parameters->seedSourceType == "multiplePremises"){
    		std::vector<std::vector<int>> sourcePremIDs;
    		// read source file in file manager (expect vector of vector of ints back)
    		read_seedSource(parameters->seedSource, sourcePremIDs);
    		std::vector<Farm*> tempPremVector;
    		for (auto& lineVector: sourcePremIDs){
    			for (auto& pId: lineVector){
    				if (farm_map.count(pId) == 0){
					  	std::cout << "Warning: Premises " << pId <<
					  	" not found, skipping this source of infection." << std::endl;
						} else {
					  	tempPremVector.emplace_back(farm_map.at(pId));
						}
    			}
    			seedFarmsByRun.emplace_back(tempPremVector);
    		}
    	}
	}
	seedFarmsByRun.swap(output);
}

/// \param[in] compare_x x-coordinate of point to be evaluated
/// \param[in] compare_y y-coordinate of point to be evaluated
/// \param[in] center_x x-coordinate of center of circle
/// \param[in] center_y y-coordinate of center of circle
/// \param[in] radius radius of circle (same units as x and y points)
/// \param[in] radiusSquared radius*radius, pre-calculated for efficiency
double Grid_manager::pointDistanceWithinRadius(const double compare_x, const double compare_y,
	const double center_x, const double center_y, const double radius, const double radiusSquared)
{
	double distanceOutput = -1; // default is not within radius (distance from center > radius)
	const double xdiff = std::abs(compare_x - center_x);
	const double ydiff = std::abs(compare_x - center_x);

	if (xdiff <= radius && ydiff <= radius){ // only bother with calculation if differences for both x and y <= radius
		double distanceSquared = xdiff*xdiff + ydiff*ydiff;
		if (distanceSquared <= radiusSquared){ // check if within circle
			distanceOutput = distanceSquared;
		}
	}
	return distanceOutput;
}

/// \param[in] compare_x x-coordinate of point to be evaluated
/// \param[in] compare_y y-coordinate of point to be evaluated
/// \param[in] center_x x-coordinate of center of circle
/// \param[in] center_y y-coordinate of center of circle
/// \param[in] radius radius of circle (same units as x and y points)
/// \param[in] radiusSquared radius*radius, pre-calculated for efficiency
unsigned int Grid_manager::pointWithinRadius(const double compare_x, const double compare_y,
	const double center_x, const double center_y, const double radius, const double radiusSquared)
{
	unsigned int output = 0;
	double pointDist = pointDistanceWithinRadius(compare_x, compare_y, center_x, center_y,
																							 radius, radiusSquared);
	if (pointDist >= 0){
		output = 1;
	}
	return output;
}

/// \param[in] cell Grid_cell to be evaluated
/// \param[in] center_x x-coordinate of center of circle
/// \param[in] center_y y-coordinate of center of circle
/// \param[in] radius radius of circle (same units as x and y points)
/// \param[in] radiusSquared radius*radius, pre-calculated for efficiency
unsigned int Grid_manager::count_cellCornersWithinRadius(Grid_cell* cell, const double center_x,
	const double center_y, const double radius, const double radiusSquared)
{
	// check if each corner is in radius
	double north = cell->Grid_cell::get_north();
	double south = cell->Grid_cell::get_south();
	double east = cell->Grid_cell::get_east();
	double west = cell->Grid_cell::get_west();

	unsigned int swInRange = pointWithinRadius(west, south, center_x, center_y, radius, radiusSquared);
	unsigned int nwInRange = pointWithinRadius(west, north, center_x, center_y, radius, radiusSquared);
	unsigned int neInRange = pointWithinRadius(east, north, center_x, center_y, radius, radiusSquared);
	unsigned int seInRange = pointWithinRadius(east, south, center_x, center_y, radius, radiusSquared);

	// return total number of corners in range
	unsigned int output = swInRange+nwInRange+neInRange+seInRange;

	if (output <= 4){
		return output;
	} else {
		std::cout << "ERROR: Grid_manager::cellCornersWithinRadius: wrong number of corners ("<< output
		<< "). Exiting..." << std::endl;
		exit(EXIT_FAILURE);
	}
}

/// Checks if neighbors in a given radius have already been determined. If not, calls
/// calc_neighborsInRadius, and returns neighbors in radius.
void Grid_manager::get_neighborsInRadius(Farm* focal, const double radius,
	const double radiusSquared, const bool distanceRequired, std::multimap<double, Farm*>& output)
{
	// check if focal Farm already has neighbors in designated radius
	double checkedRadius = focal->Farm::get_neighborRadiusCalculated();

	if (checkedRadius < radius){
	// need to generate multimap of neighbors/distances
		std::multimap<double, Farm*> dn; // blank multimap to be written to
		calc_neighborsInRadius(focal, radius, radiusSquared, distanceRequired, dn);
		dn.swap(output);
	} else {
	// neighbors have already been calculated at this radius or larger
		// if distancesRequired, check for missing distance values (may happen when all prems in a cell are added)
		if (distanceRequired == 1){
			const std::multimap<double, Farm*>* dnCheck = focal->Farm::get_distancesNeighbors();
			auto missingDistance = dnCheck->equal_range(-1);
			// if any neighbors are missing distances, fill in
			if (std::distance(missingDistance.first, missingDistance.second) > 0){
				fill_premisesDistances(focal);
			}
		}
		// point to focal premises' distances-neighbors
		const std::multimap<double, Farm*>* dn = focal->Farm::get_distancesNeighbors();

		// if neighbors calculated from larger radius, filter to this radius value
		if (checkedRadius > radius){
			auto endRange = dn->upper_bound(radius);
			std::multimap<double, Farm*> tempOutput (dn->begin(), endRange);
			tempOutput.swap(output);
		} else if (checkedRadius == radius){
			std::multimap<double, Farm*> tempOutput (*dn); // copy as-is
			tempOutput.swap(output);
		}
	}

}

/// \param[in] focal Focal premises
/// \param[in] radius Radius (in same units as premises coordinates) to use as threshold
/// \param[in] radiusSquared radius*radius, pre-calculated for efficiency
/// \param[out] output Vector to which premises within radius will be copied
void Grid_manager::calc_neighborsInRadius(Farm* focal, const double radius,
	const double radiusSquared, bool distanceRequired, std::multimap<double, Farm*>& output)
{
	const int focalCellID = focal->Farm::get_cellID();
	const double center_x = focal->Farm::get_x();
	const double center_y = focal->Farm::get_y();

	std::vector<Grid_cell*> cellsToCheck;
	cellsToCheck.emplace_back(allCells.at(focalCellID));
	int index = 0;

	std::multimap<double, Farm*> distancesNeighbors;

	while (index < cellsToCheck.size()){
		unsigned int numCornersInRadius = count_cellCornersWithinRadius(cellsToCheck.at(index),
			center_x, center_y, radius, radiusSquared);
		std::vector<Farm*> inCell = cellsToCheck.at(index)->Grid_cell::get_farms();

		// if all corners are in radius, add all farms in cell as neighbors
		if (numCornersInRadius ==4){
			if (distanceRequired == 0){
				add_premisesWithoutDistance(inCell, distancesNeighbors);
			} else {
				add_premisesInRadius(center_x, center_y, radius, radiusSquared, inCell, distancesNeighbors);
			}
			// get neighbors of this cell and add to queue
			get_neighborCellsToCheck(cellsToCheck.at(index), cellsToCheck);
		// if there are 1-3 corners in radius, check each farm in cell
		} else if (numCornersInRadius >0 && numCornersInRadius <4){
			// check if each farm in cell is in radius
			add_premisesInRadius(center_x, center_y, radius, radiusSquared, inCell, distancesNeighbors);
			// get neighbors of this cell and add to queue
			get_neighborCellsToCheck(cellsToCheck.at(index), cellsToCheck);
		// if there are no corners in radius, but this is the focal cell, check all farms in cell
		} else if (numCornersInRadius == 0 && cellsToCheck.at(index)->Grid_cell::get_id() == focalCellID){
			// cell corners are outside of radius but since this is the focal cell, check all farms within
			add_premisesInRadius(center_x, center_y, radius, radiusSquared, inCell, distancesNeighbors);
			// get neighbors of this cell and add to queue
			//get_neighborCellsToCheck(cellsToCheck.at(index), cellsToCheck);
		}
		// otherwise, if numCornersInRadius is 0 and not focal cell, do nothing
		index++;
	}

	output.swap(distancesNeighbors);

	// can comment out next two lines for less memory usage
	focal->Farm::set_distancesNeighbors(distancesNeighbors);
	focal->Farm::set_neighborRadiusCalculated(radius);
}

/// \param[in] focalCell Grid_cell for which to retrieve neighboring cells
/// \param[in] cellsToCheck Vector of cells to be checked, to which neighboring cells will be copied
void Grid_manager::get_neighborCellsToCheck(Grid_cell* focalCell, std::vector<Grid_cell*>& cellsToCheck)
{
	const std::vector<Grid_cell*>* neighborCells = focalCell->Grid_cell::get_neighbors();
	for (auto& nc:(*neighborCells)){
		auto checkIfAdded = std::find(cellsToCheck.begin(), cellsToCheck.end(), nc);
		if (checkIfAdded == cellsToCheck.end()){ // is not already in list
			cellsToCheck.emplace_back(nc);
		}
	}
}

/// \param[in] center_x x-coordinate of center of circle
/// \param[in] center_y y-coordinate of center of circle
/// \param[in] radius radius of circle (same units as x and y points)
/// \param[in] radiusSquared radius*radius, pre-calculated for efficiency
/// \param[in] premisesToCheck Vector of premises to check for being within radius
/// \param[in] neighbors Vector to which premises that are within radius are copied
void Grid_manager::add_premisesInRadius(const double center_x, const double center_y,
	const double radius, const double radiusSquared, const std::vector<Farm*>& premisesToCheck,
	std::multimap<double, Farm*>& distancesNeighbors)
{
	for (auto& f:premisesToCheck){
		double premisesDistance = pointDistanceWithinRadius(f->Farm::get_x(), f->Farm::get_y(),
																										center_x, center_y, radius,
																										radiusSquared);
		if (premisesDistance >= 0){
			distancesNeighbors.emplace(premisesDistance, f);
		}
	}
}

/// Adds premises as a neighbor with temporary distance of -1 (calculated later if needed)
void Grid_manager::add_premisesWithoutDistance(const std::vector<Farm*>& premisesToAdd,
	std::multimap<double, Farm*>& distancesNeighbors)
{
	for (auto& f:premisesToAdd){
		distancesNeighbors.emplace(-1,f);
	}
}


/// Returns neighboring cells containing premises in the same state. If
/// inAnyReportedStates is TRUE, also returns any neighboring cells from states in which
/// disease has been reported.
/// \param[in] focalCell Grid_cell for which to retrieve neighboring cells
/// \param[in] focalState State containing premises of interest
/// \param[in] reportedStates Vector of (2 letter abbreviations of) reported states
/// \param[in] inAnyReportedStates Boolean value indicating if cells in other states that have been reported should be included as neighboring cells
void Grid_manager::get_neighborCellsByState(Grid_cell* focalCell, std::string focalState,
	std::vector<std::string>& reportedStates, bool inAnyReportedStates,
	std::vector<Grid_cell*>& output)
{
	const std::vector<Grid_cell*>* neighborCells = focalCell->Grid_cell::get_neighbors();
	for (auto& nc:(*neighborCells)){
		const std::set<std::string> statesInNeighborCell = nc->Grid_cell::get_states();
		auto /*pointer*/ cellHasFocalState = std::find(statesInNeighborCell.begin(), statesInNeighborCell.end(), focalState);
		if (cellHasFocalState != statesInNeighborCell.end()){ // if neighbor cell contains focal state
			output.emplace_back(nc); // remove duplicates when comparing to other cells' neighbors
		} else if (inAnyReportedStates==1){ // cell does not contain focal state
			// check if cell includes a reported state
			auto inReported = reportedStates.end();
			auto s = statesInNeighborCell.begin();
			while (inReported == reportedStates.end() && s != statesInNeighborCell.end()){
				inReported = std::find(reportedStates.begin(), reportedStates.end(), *s);
				s++;
			}
			if (inReported != reportedStates.end()){ // cell contains reported state, add to list
				output.emplace_back(nc);
			}
		}
	}
}



/// This function should only be called when distances to each premises are needed (i.e.
/// to prioritize by proximity) and there are unknown distances to fill in, for efficiency
void Grid_manager::fill_premisesDistances(Farm* focal)
{
	double focal_x = focal->Farm::get_x();
	double focal_y = focal->Farm::get_y();
	// make a copy of the distances-neighbors multimap (includes some farms without distances)
	std::multimap<double, Farm*> dn = *(focal->Farm::get_distancesNeighbors());
	// get range of elements where distance == -1
	auto itRange = dn.equal_range(-1);

	for (std::multimap<double, Farm*>::iterator it = itRange.first; it!=itRange.second; ++it){
		Farm* neighborPrem = it->second;
		double neighbor_x = neighborPrem->Farm::get_x();
		double neighbor_y = neighborPrem->Farm::get_y();
		double xdiff = focal_x - neighbor_x;
		double ydiff = focal_y - neighbor_y;
		double distSquared = xdiff*xdiff + ydiff*ydiff;
		dn.emplace(distSquared, neighborPrem); // add new distance calculation
	}
	// remove previous versions of elements without distance
	dn.erase(itRange.first, itRange.second);
	// assign modified multimap to Farm object
	focal->Farm::set_distancesNeighbors(dn);
}


/// Finds cell containing a point (i.e. landfill) with x and y coordinates and FIPS code.
/// If x, y, and FIPS code don't align, will return -1. If FIPS code not defined from
/// premises list, will return -1. Could potentially add code in the future to assign
/// cells in these cases. Called when Control resources are created to determine which
/// cell it's within.
/// \param[in] x X-coordinate of point
/// \param[in] y Y-coordinate of point
/// \param[in] countyID FIPS code of point
int Grid_manager::get_parentCell(double x, double y, std::string countyID){
	if (cellsByCounty.count(countyID) > 0){
		std::vector<Grid_cell*> possibleCells = cellsByCounty.at(countyID);
		if (possibleCells.size() == 1){
			return possibleCells.at(0)->Grid_cell::get_id();
		} else {
			for (auto& cell:possibleCells){
				// if x/y are in cell, return that cell id
				if (x <= cell->Grid_cell::get_east() && x >= cell->Grid_cell::get_west() &&
						y <= cell->Grid_cell::get_north() && y >= cell->Grid_cell::get_south()){
					return cell->Grid_cell::get_id();
				}
			} // if no matches found, will get past this for-statement; return -1
			return -1;
		}
	} else {
		return -1;
	}
}

// returns the norminf value for the species
const std::unordered_map<std::string, double>& Grid_manager::get_normInf_map() const{

    return normInf;
}

const std::unordered_map<std::string, double>& Grid_manager::get_normSus_map() const{

    return normSus;
}

