#ifndef Shipment_manager_h
#define Shipment_manager_h

#include <set>
#include <tuple> // std::tuple

#include "Grid_manager.h"
#include "shared_functions.h"
#include <gsl_rng.h>

extern int verboseLevel;

class County;
class Status_manager;

struct Shipment // used in Shipment, Status
{
	int timestep; ///< Timestep of shipment
	size_t day_of_year; ///< Day of the year of the shipment.
	int origID; ///< Premises ID of shipment origin
	int destID; ///< Premises ID of shipment destination
	std::string origFIPS; ///< County ID of shipment origin
	std::string destFIPS; ///< County ID of shipment destination
	std::string species; ///< Species or animal type in shipment
	size_t volume; ///< Number of animals shipped.
	bool infectious; ///< T/F: this is from an infectious premises
	int ban; ///< Ban prevented shipment? 0 = no, 1 = yes
	std::string time_period; ///Time period of shipment. From config file option 45.
	int origState_id;
	std::string origState_abbrev;
	int destState_id;
	std::string destState_abbrev;
};

struct coShipment{
	int t; ///< Time of shipment
	std::string origFIPS; ///< Premises ID of shipment origin
	std::string destFIPS; ///< Premises ID of shipment destination
	std::string species;
	int volume; ///> Generally 1, but allows for multiple shipments between counties in the same t
	std::string ban; ///< Ban level: '' = no ban, 'implemented', 'effective', 'inactive'
};

///> Manages the shipments (USAMM) part of the simulation
/// Predicts shipments from counties with infectious farms
/// Gets called at each timepoint
/// Bans: go into effect "after" shipments determined - for tracking business continuity
class Shipment_manager
{
	private:
		int verbose; ///< Can be set to override global setting for console output
        bool shipments_on;
		// const pointers to Grid_manager objects, parameters:
		const std::unordered_map<std::string, County*>* FIPS_map;
		const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpeciesMap;
		const Parameters* parameters;

		Status_manager* S; ///< Const pointer to Status manager (to access up-to-date premises statuses)
		// used to generate random shipments
		std::vector<County*> allCounties; // list of all possible destination FIPS, based on premises file
		std::set<State*> allStates_set; //Set containing all states.
		int farmFarmMethod;
		std::vector<std::string> species;
		std::unordered_map<std::string, std::vector<std::string>> speciesFIPS; // just for countycounty fake assignment to make sure appropriate county is chosen

        gsl_rng* R; //A gsl random number generator. Initialized in initialize()
		// the following are recreated/rewritten at each timestep
		std::vector<coShipment>
			countyShipmentList; // coShipment defined above
		std::vector<Shipment*> // shipment defined in shared_functions.h
			farmShipmentList;
		int startRecentShips, startCoRecentShips; // indicates index in shipmentList where the most recent set of shipments starts

		// functions
		void initialize();
		Farm* largestStatus(std::vector<Farm*>&, std::string&); ///< Finds largest premises with "status", from vector sorted by population

		///Creates and returns a pointer to a shipment struct.
		Shipment* generateInfectiousShipment(Farm* origin_farm, size_t timestep, size_t day_of_year,
                                             const std::string& time_period);

	public:
		Shipment_manager(
			const std::unordered_map<std::string, County*>* in_FIPS_map, // a map of FIPS codes to farms
			const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpMap, // sorted populations of species on farms
			Status_manager* in_S,
			int ffm, // farm assignment method
			const std::vector<std::string>& speciesOnPrems, // list of species on premises
			const Parameters* p);

        //Constructor without status manager. Used for network generation since that does not require keeping track of statuses.
        Shipment_manager(
			const std::unordered_map<std::string, County*>* in_FIPS_map, // a map of FIPS codes to farms
			const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpMap, // sorted populations of species on farms
			int ffm, // farm assignment method
			const std::vector<std::string>& speciesOnPrems, // list of species on premises
			const Parameters* p);

		~Shipment_manager();

        ///Randomly creates shipments originating from infected farms during one time step. If the
        ///infFarm vector argument is empty will create a complete shipment network for a whole year.
		void makeShipmentsMultinomial(size_t timestep, size_t days_in_period, size_t days_rem,
                                      std::string time_period, std::vector<Shipment*>& output,
                                      std::vector<Farm*>& infFarms, std::vector<Farm_type*> ft_vec);
		///Generates and writes to file a complete yearly network of shipments based on the shipment parameters.
		void makeNetwork(std::vector<std::string> out_fname, Grid_manager& G);

		std::string formatOutput(int, int); // formats output to string

};

#endif
