#ifndef Grid_manager_h
#define Grid_manager_h

#include "File_manager.h" // for parameter struct
#include "Grid_cell.h"
#include "shared_functions.h" //random_unique
#include "USAMM_parameters.h"

#include <algorithm> // std::sort, std::any_of, std::find
#include <map> // std::multimap
#include <stack>
#include <tuple>
#include <unordered_map>
#include <utility> // std::pair


class County;
class State;

extern int verboseLevel;

///  Creates set of Grid_cells determined by local farm density, stores relevant values with Farms
class Grid_manager
{
	private:
		int verbose; ///< Can be set to override global setting for console output

		const Parameters* parameters;
		// variables for grid creation
		unsigned int maxFarms; ///< Threshold number of premises per cell (cell size takes precedence)
		bool shipments_on; ///Keeps track of if shipments are turned on.

		std::string shipment_kernel_str; ///Stores the shipment kernel type.
		std::unordered_map<int, Grid_cell*>
			allCells; ///< Unordered_map of all cells in grid
		std::unordered_map<int, Farm*>
			farm_map; ///< Unordered_map of all premises objects
        std::vector<Farm*> farm_vector; ///<Vector containing pointers to all premises objects.
		std::unordered_map<std::string, State*>
            state_map; //Contains states, name as key.
        std::vector<State*> state_vector;
		std::unordered_map<std::string, County*>
			FIPS_map; // key is fips code, value is county object
    std::vector<County*>
            FIPS_vector;
		std::unordered_map<std::string,
        std::unordered_map< std::string, std::vector<Farm*> >> fipsSpeciesMap;
			// key is fips code, then species name, then sorted by population size
		std::unordered_map<std::string, std::vector<Grid_cell*>> cellsByCounty;
 		std::vector<Farm*>
 			farmList; // vector of pointers to all farms (deleted in chunks as grid is created)

		std::tuple<double,double,double,double>
			xylimits; ///< Ranges of premises coordinates: [0]x min, [1]x max, [2]y min, [3]y max

		// variables for infection evaluation
		std::vector<std::string> speciesOnPrems; ///< List of species on all farms provided in premises file
		std::unordered_map<std::string,double> susExponents; ///< Species-specific susceptibility exponents, in same order as speciesOnAllFarms
		std::unordered_map<std::string,double> infExponents; ///< Species-specific infectiousness exponents, in same order as speciesOnAllFarms
		std::unordered_map<std::string,double> susValues; ///< Species-specific susceptibility values, in same order as speciesOnAllFarms
		std::unordered_map<std::string,double> infValues; ///< Species-specific infectiousness values, in same order as speciesOnAllFarms
		std::unordered_map<std::string,double> normInf; ///< Normalized species-specific infectiousness values, in same order as speciesOnAllFarms
		std::unordered_map<std::string,double> normSus; ///< Normalized species-specific susceptibility values, in same order as speciesOnAllFarms
		Local_spread* kernel;

		unsigned int committedFarms; ///< Used to double-check that all loaded premises were committed to a cell
		int printCellFile;
		std::string batch; ///< Cells printed to file with name: [batch]_cells.txt
		int USAMM_temporal_index;
		int start_day;
		size_t days_rem_of_period = 0;
		size_t days_in_period = 0;
		size_t USAMM_current_year = 1;
		std::string USAMM_temporal_name;
		std::map<Farm_type*, USAMM_parameters> usamm_parameters;
		std::unordered_map<std::string, Farm_type*> farm_types_by_herd;
		std::unordered_map<std::string, Farm_type*> farm_types_by_name;
		std::vector<Farm_type*> farm_types_vec;

		// functions
		///Reads counties and states from file specified in config #18.
		void readFips_and_states();
		///Does all the reading of premises-file (FLAPS) specified in config #11
		void readFarms(const std::string& farm_fname);
		///Sets the covariates of the counties.
		void initFipsCovariatesAndCounties();
		///Calculates the county-level shipping probabilities based on the
		///current state-level parameters and
		///nation-level covariate parameters given the time period.
		void updateFipsShipping(std::string time_period);
		///Updates the shipment rate for the states using new farm and county covariate weights
		///as well as the latest lambda_null.
		void updateStateLambdas();
		///Update all the states with new shipping parameters (std, kurtosis, N and s)
		///for the given time period.
		void updateStatesShipping(std::string time_period, size_t days_in_period, size_t days_rem, bool reset);
		void set_maxFarms(unsigned int in_maxFarms); //inlined
		std::string to_string(Grid_cell&) const;
		std::vector<Farm*> getFarms(
			std::tuple<int,double,double,double>& cellSpecs,
			const unsigned int maxFarms=0); ///< Makes list of farms in a cell (quits early if over max)
		void removeParent(
			std::stack< std::tuple<int,double,double,double> >& queue);///< Removes 1st item in queue
		void addOffspring(
			std::tuple<int,double,double,double> cellSpecs,
			std::stack< std::tuple<int,double,double,double> >& queue);
		void commitCell(
			std::tuple<int,double,double,double> cellSpecs,
			std::vector<Farm*>& farmsInCell); ///< Adds Grid_cell to allCells
		void splitCell(
			std::tuple<int,double,double,double>& cellSpecs,
			std::stack< std::tuple<int,double,double,double> >& queue); ///< Replaces parent cell with subdivided offspring quadrants
 		void assignCellIDtoFarms(int cellID, std::vector<Farm*>& farmsInCell);
 		void removeFarmSubset(std::vector<Farm*>&, std::vector<Farm*>&); ///< Remove farms in first vector from second vector

		void read_seedSource(std::string, std::vector<std::string>&); ///< Reads seed file with one FIPS code per line
		void read_seedSource(std::string, std::vector<int>&); ///< Reads seed file with one premises ID per line
		void read_seedSource(std::string, std::vector<std::vector<int>>&); ///< Reads seed file with multiple premises IDs per line
		void select_randomPremisesPerCounty(std::vector<std::vector<Farm*>>&); ///< Selects seed premises from all counties
		void select_randomPremisesPerCounty(std::vector<std::string>, std::vector<std::vector<Farm*>>&); ///< Selects seed premises from specified counties
		double pointDistanceWithinRadius(const double, const double, const double, const double, const double, const double); ///< Returns distance of x and y coordinates from center point IF within radius (if not, returns -1)
		unsigned int pointWithinRadius(const double, const double, const double, const double, const double, const double); ///< Boolean wrapper function for pointDistanceWithinRadius (used for checking cell corners, where distance need not be recorded), returns 1 if within radius, 0 if not
		unsigned int count_cellCornersWithinRadius(Grid_cell*, const double, const double, const double, const double); ///< Evaluates number of corners of grid cell within radius, returns 0-4
		void calc_neighborsInRadius(Farm*, const double, const double, bool, std::multimap<double, Farm*>&); ///< Retrieves premises within radius from focal premises
		void get_neighborCellsToCheck(Grid_cell*, std::vector<Grid_cell*>&); ///< Adds neighboring cells to vector of cells to check
		void add_premisesInRadius(const double, const double, const double, const double,
			const std::vector<Farm*>&, std::multimap<double, Farm*>&); ///< From vector of premises to check, adds premises that are within radius to multimap of neighbors
		void add_premisesWithoutDistance(const std::vector<Farm*>&,
			std::multimap<double, Farm*>&); ///< Adds premises as neighbors without specifying distance
		void fill_premisesDistances(Farm*); ///< Fill in missing distances to other premises

		// functions for infection evaluation
		double shortestCellDist2(Grid_cell*, Grid_cell*); ///< Calculates (shortest distance between two cells)^2
		void makeCellRefs(); ///< Calculates and stores kernel values and other pre-processing tasks
		// functions for infection evaluation
		void set_FarmSus(Farm*); ///< Calculates premises susceptibility and stores in Farm
		void set_FarmInf(Farm*); ///< Calculates premises infectiousness and stores in Farm


	public:
		Grid_manager(const Parameters*);
		~Grid_manager();

		// main function that splits/commits cells:
		// 1st way to initiate a grid: specify maximum farms per cell, min size
		void initiateGrid(
			const unsigned int,
			const int);

		// 2nd way to initiate a grid: specify file containing cell specs
		void initiateGrid(
			std::string &cname);

		// 3rd way to initiate a grid: specify side length (same units as x/y) for uniform cells
		void initiateGrid(
			double cellSide);

		const std::unordered_map<int, Grid_cell*>*
			get_allCells() const; //inlined

		const std::unordered_map<int, Farm*>*
			get_allFarms() const; //inlined

        std::vector<Farm*>& get_allFarms_vector(); //Inlined

		const std::unordered_map<std::string, County*>*
			get_allCounties() const; //inlined

		const std::unordered_map<std::string, State*>*
			get_allStates() const; //inlined

		const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>*
			get_fipsSpeciesMap() const; //inlined

		void get_seedPremises(std::vector<std::vector<Farm*>>&);

		void printCells();
		void get_neighborsInRadius(Farm*, const double, const double, const bool,
			std::multimap<double, Farm*>&); ///< Return neighboring premises within radius
		void get_neighborCellsByState(Grid_cell*, std::string,
			std::vector<std::string>&, bool, std::vector<Grid_cell*>&);
		int get_parentCell(double, double, std::string);

        ///Constructs the USAMM_parameters objects and assigns shipping parameters
        ///to states and covariates to counties.
        void initShippingParameters(int t, int start_day_in);
        ///Checks to see if a transition has been made from one time period to the next
        ///since the last time step (eg. Q1 -> Q2) and updates all the USAMM parameters
        ///if necessary. This function must be called each timestep of the simulation.
		void updateShippingParameters(int t, bool restart = false);
		///Interface to get the USAMM .res file line that the farm-type
		///specific USAMM_parameter object has loaded.
		std::string get_generation_string(Farm_type* ft);

        ///Sums the origin shipping weight for all counties within a state and normalizes
        ///each countys weight with this norm so that they sum to 1 within a given state.
        ///Does this for each state.
		void normalizeShippingWeights();
		///Return the current time period that the simulatin is currently in.
		std::string get_time_period();
		///Returns the number of days in the current period.
		size_t get_days_in_period();
		///Returns the remaining number of days in the current period.
		size_t get_rem_days_of_period();

        Farm_type* get_farm_type_by_herd(std::string herd); //Returns pointer to type based on what species are present on a farm.
		Farm_type* get_farm_type_by_name(std::string name);
		std::vector<Farm_type*> get_farm_types();

        const std::unordered_map<std::string, double>& get_normInf_map() const;
        const std::unordered_map<std::string, double>& get_normSus_map() const;
};

/// "Compare" function to sort farms by x-coordinate,
/// used when assigning farms to uniform cell grid
inline bool sortByX(const Farm* farm1, const Farm* farm2)
{
	return (farm1 -> get_x()) < (farm2 -> get_x());
}

inline void Grid_manager::set_maxFarms(unsigned int in_maxFarms)
{
	maxFarms = in_maxFarms;
}

inline const std::unordered_map<int, Grid_cell*>*
	Grid_manager::get_allCells() const
{
	return &allCells;
}

inline const std::unordered_map<int, Farm*>*
	Grid_manager::get_allFarms() const
{
	return &farm_map;
}

inline std::vector<Farm*>& Grid_manager::get_allFarms_vector()
{
    return farm_vector;
}

inline const std::unordered_map<std::string, County*>*
	Grid_manager::get_allCounties() const
{
	return &FIPS_map;
}
inline const std::unordered_map<std::string, State*>*
	Grid_manager::get_allStates() const
{
	return &state_map;
}
inline const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>*
	Grid_manager::get_fipsSpeciesMap() const
{
	return &fipsSpeciesMap;
}

// used to look up re-used cell distances
template<typename T> std::vector<T> orderNumbers(T& number1, T& number2)
// order number1 and number2 from lowest to highest
{
	std::vector<T> ordered;
	ordered.emplace_back(number1);
	if (number2 < number1){
		ordered.insert(ordered.begin(),number2);
	} else {
		ordered.emplace_back(number2); // if number2 is larger or equal to number1
	}
	return ordered;
}

/// Sorts farms by population for a given species/type
struct comparePop
{
	comparePop(std::string in_species) : species(in_species) {}
	bool operator() (const Farm* farm1, const Farm* farm2){
		return (farm1->Farm::get_size(species)) < (farm2->Farm::get_size(species));
	}
	private:
		std::string species;
};


#endif
