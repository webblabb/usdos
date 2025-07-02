//© 2024 Colorado State University
#ifndef Shipment_manager_h
#define Shipment_manager_h

#include <set>
#include <tuple> // std::tuple

#include "Grid_manager.h"
#include "shared_functions.h"
#include <gsl/gsl_rng.h>

extern int verboseLevel;

class County;
class Status_manager;

typedef std::pair<County*, std::vector<Farm*>> c_farmvec_pair;
typedef std::vector<c_farmvec_pair> Vec_cf_1d;
typedef std::vector<Vec_cf_1d> Vec_cf_2d;
typedef std::vector<Vec_cf_2d> Vec_cf_3d;

typedef std::vector<Farm*> Vec_fp_1d;
typedef std::vector<Vec_fp_1d> Vec_fp_2d;
typedef std::vector<Vec_fp_2d> Vec_fp_3d;
typedef std::vector<Vec_fp_3d> Vec_fp_4d;

struct Shipment // used in Shipment, Status
{
	int timestep; ///< Timestep of shipment
	size_t day_of_year; ///< Day of the year of the shipment.
	Farm* oPrem; ///< Pointer to shipment origin farm.
	Farm* dPrem; ///< Pointer to shipment destination farm.
	std::string origFIPS; ///< County ID of shipment origin
	std::string destFIPS; ///< County ID of shipment destination
	std::string species; ///< Species or animal type in shipment
	int volume; ///< Number of animals shipped.
	int n_infected; ///< The number of infected animals on the shipment.
	int ban; ///< Ban prevented shipment? 0 = no, 1 = yes
	std::string time_period; ///Time period of shipment. From config file option 45.
	int origState_id;
	std::string origState_abbrev;
	int destState_id;
	std::string destState_abbrev;
	std::string originIndType; ///< Industry type of origin premises (farm, feedlot, market).
	int originSize; ///< Size of origin premises.
	std::string destIndType; ///< Industry type of destination premises (farm, feedlot, market).
	int destSize; ///< Size of destination premises.
	Prem_status* oPst; ///<The corresponding prem status of the origin farm if there is one.
	std::string permitType;
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
	    std::vector<County*> allCounties; // list of all possible destination FIPS, based on premises file
        const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpeciesMap;
		const Parameters* parameters;
		std::vector<std::string> species;

		int verbose; ///< Can be set to override global setting for console output
        bool shipments_on;
		// const pointers to Grid_manager objects, parameters:
		const std::unordered_map<std::string, County*>* FIPS_map;
		// used to generate random shipments
		std::set<State*> allStates_set; //Set containing all states.
		std::unordered_map<std::string, std::vector<std::string>> speciesFIPS; // just for countycounty fake assignment to make sure appropriate county is chosen

        gsl_rng* R; //A gsl random number generator. Initialized in initialize()
		// the following are recreated/rewritten at each timestep
		std::vector<coShipment>
			countyShipmentList; // coShipment defined above
		int startRecentShips, startCoRecentShips; // indicates index in shipmentList where the most recent set of shipments starts
        Vec_d_2d p_obs_COM_lookup_table;
		// functions
		void initialize();

		///Creates and returns a pointer to a shipment struct.
		Shipment* generateShipmentUSAMMv2(Farm* origin_farm, size_t timestep, size_t day_of_year,
                                          const std::string& time_period);

        size_t make_cattle_shipment_volume(Farm* o_prem, Farm* d_prem, Prem_class* o_pcl, Prem_class* d_pcl,
                                           std::string time_period, USAMMv3_parameters& up, Prem_class* mkt_pcl);
        size_t make_swine_shipment_volume(Farm* o_prem, Farm* d_prem, Prem_class* o_pcl, Prem_class* d_pcl,
                                          std::string time_period, USAMMv3_parameters& up);

	public:
		Shipment_manager(
			const std::vector<County*> in_FIPS_vec, // a map of FIPS codes to farms
			const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpMap, // sorted populations of species on farms
			const Parameters* p);

		~Shipment_manager();

        ///Randomly creates shipments originating from infected farms during one time step. If the
        ///infFarm vector argument is empty will create a complete shipment network for a whole year.
		void makeShipmentsUSAMMv2(size_t timestep, size_t days_rem,
                                      std::string time_period, std::vector<Shipment*>& output,
                                      std::vector<Farm*>& infFarms, std::vector<Farm_type*> ft_vec);

        void makeShipmentsUSAMMv3(size_t timestep, size_t day_of_year, std::string time_period,
                                  int time_period_idx, std::vector<Shipment*>& output,
                                  std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec);
        void makeShipmentsUSAMMv3_new(size_t timestep, size_t day_of_year, std::string time_period,
                                  int time_period_idx, std::vector<Shipment*>& output,
                                  std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec);

        void makeNetworkUSAMMv3(size_t timestep, size_t day_of_year, std::string time_period,
                                int time_period_idx, std::vector<Shipment*>& output,
                                std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec);
		///Generates and writes to file a complete yearly network of shipments based on the shipment parameters.
		void makeNetworkUSAMMv2(std::vector<std::string> out_fnames, Grid_manager& G);

		void makeNetworkUSAMMv3(std::vector<std::string> out_fnames, Grid_manager& G);

		std::string formatOutput(int, int); // formats output to string

};

#endif
