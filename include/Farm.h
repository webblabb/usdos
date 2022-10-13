#ifndef FARM_H
#define FARM_H

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <map> // for distances to neighbors
#include <string>
#include <unordered_map>
#include <utility> // for std::iter_swap in rem_probPreventExposure
#include <vector>
#include "Point.h"

class County;
class State;
class Farm_type;
class Shipment;
struct BtbWithinHerdParams;
struct Parameters;
struct Prem_class;

extern int verboseLevel;

/// Describes a premises - one of these objects is created for each premises row read in
/// from the premises file.
class Farm
{
	protected: // allows access from derived class Prem_status
		int id,	///< Unique integer identifier read from premises file
			cellID; ///< Integer identifier of grid_cell assigned to this premises during grid creation
		double x_coordinate, ///< x-coordinate from projected longitude (same units as local spread kernel)
			y_coordinate, ///< y-coordinate from projected latitude (same units as local spread kernel)
			inf_max, ///< The maximum attainable infectiousness for this premises, same as inf if partial transition is off, peak value if partial transition is on.
			sus_max; ///Highest susceptibility over the year when using quarterly premises sizes.
        std::vector<double> sus, ///< Calculated total susceptibility of this premises by quarter
                            inf; ///< Calculated total infectiousness of this premises by quarter
        Point position;
		County* parent_county;
		int idx_in_county; //This prem is stored at this index of the vector of Farm object pointers in the parent county. Used in the shipment generation code.
		State* parent_state;
		Farm_type* farm_type;
		std::string farm_type_str;
		std::string fips; ///< County identifier (FIPS code)
		const Parameters* parameters;
		const BtbWithinHerdParams* btb_p;
		bool is_frm = false;
		bool is_fdl = false;
        bool is_mkt = false;
        double slaughter_shipment_factor = -1.0; ///< Relates this premises' rate of shipping to slaughter to its "regular" shipment.
		double yearly_slaughter_fraction = -1.0; ///< The proportion of this premises' inventory that goes to slaughter every year. Can be larger than 1.0.
		double latest_shipping_rate = -1.0; ///< The total shipment rate from this premises according to last update of USAMMv3 parameters. Only set if shipments have been genereated from this premises.
		Prem_class* prem_class;
		double neighborRadiusCalculated; ///< If neighbors have been calculated and stored as distancesNeighbors, the radius searched (else 0)
//		std::unordered_map< std::string, int > speciesCounts; ///< Numbers of animals of each type currently in use, keyed by types
		std::vector< std::vector<int> > speciesCountsQuarterly; ///<Number of animals by farm type index and quarter.
		std::vector< std::unordered_map<std::string, int>> speciesCountsMapQuarterly; ///<Number of animals by quarter and species string.
		std::multimap<double, Farm*> distancesNeighbors; ///< Stored calculations of distances to neighboring premises
        int USAMMv3_size = 0; ///< The total number of animals on the premises. Used in the USAMMv3 shipment generation process.
        int USAMMv3_binned_size = 0; ///< The total number of animals on the premises binned according to the USAMMv3 premises size binning function.
		int USAMMv3_unbinned_size_idx = 2;
		std::map<Farm_type*, double> oweight; ///< The origin weight of this farm in relation to all other farms within the same state of the same type.
		std::map<Farm_type*, double> normalized_oweight; ///< The origin weight of this farm in relation to all other farms within the same state of the same type, normalized so the sum of all farms in state = 1.0.
		std::map<Farm_type*, double> dweight; ///< The destination weight of this farm in relation to all other farms within the same state of the same type.
		std::map<Farm_type*, double> normalized_dweight; ///< The destination weight of this farm in relation to all other farms within the same state of the same type, normalized so the sum of all farms in state = 1.0.

		static size_t current_quarter_idx;

	public:
		Farm(int, double, double, std::string, const Parameters*);
		~Farm();
		int get_id() const; // inlined
		int get_cellID() const; // inlined
		Farm_type* get_farm_type() const; //inlined
		Prem_class* get_prem_class() const { return prem_class; }
		int get_prem_class_idx() const;
		double get_x() const; // inlined
		double get_y() const; // inlined
		double get_sus(unsigned int quarter_idx) const; // inlined
		double get_inf(unsigned int quarter_idx) const; // inlined
        double get_sus_max() const; ///< Returns maximum susceptibility across quarters.
        double get_inf_max() const; ///< Returns maximum infectiousness where the partial transmission function peaks across quarters.
		int get_USAMMv3_binned_size() const { return USAMMv3_binned_size; }
		int get_USAMMv3_unbinned_size() const;
        double get_normalized_oweight(Farm_type* ft) const;
		double get_unnormalized_oweight(Farm_type* ft) const;
		double get_unnormalized_dweight(Farm_type* ft) const;
		double get_shipment_rate();
		double get_slaughter_shipment_factor() { return slaughter_shipment_factor; }
		double get_yearly_slaughter_fraction() { return yearly_slaughter_fraction; }
 		std::string get_fips() const; // inlined

 		bool is_farm() const { return is_frm; };
 		bool is_feedlot() const { return is_fdl; }
 		bool is_market() const { return is_mkt; }
 		const std::unordered_map<std::string, int>& get_spCountsCurrentQuarter(); // inlined
 		const std::unordered_map<std::string, int>& get_spCountsSpecificQuarter(int quarter_idx); // inlined
 		int get_size_current_quarter(const Farm_type* ft) const;
 		int get_size_specific_quarter(const Farm_type* ft, int quarter_idx) const;
 		int get_size_allSpecies() const;
		County* get_parent_county() const; //Inlined
		int get_idx_in_county() { return idx_in_county; }
		void set_idx_in_county(int idx) { idx_in_county = idx; }
 		State* get_parent_state() const; //Inlined
 		const std::multimap<double, Farm*>* get_distancesNeighbors(); //inlined
 		double get_neighborRadiusCalculated() const; //inlined
 		void set_xy(double x, double y);
		void set_cellID(const int cellID);
		void set_farm_type(Farm_type* in_type);
 		void set_quarterlySpeciesCounts(int ft_idx, std::string ft_name, std::vector<int>); //Sets the initial vector of sizes for each quarter.
		void set_slaughter_shipment_factor(double factor) { slaughter_shipment_factor = factor; }
		void set_yearly_slaughter_fraction(double fraction) { yearly_slaughter_fraction = fraction; }
		void set_prem_class(Prem_class* in_class);
 		void set_sus(const double sus, unsigned int quarter_idx);
 		void set_sus_all_quarters(const double);
 		void set_inf(const double inf, unsigned int quarter_idx);
 		void set_inf_all_quarters(const double);
 		void set_USAMMv3_size(int val) { USAMMv3_size = val; }
 		void set_USAMMv3_binned_size(int val) { USAMMv3_binned_size = val; }
 		void set_inf_max(const double in_inf_max);
 		void set_normalized_oweight(const double in_oweight, Farm_type* ft);
 		void set_normalized_dweight(const double in_dweight, Farm_type* ft);
 		void set_unnormalized_oweight(const double in_oweight, Farm_type* ft);
 		void set_unnormalized_dweight(const double in_dweight, Farm_type* ft);
 		void set_latest_shipping_rate(double rate) { latest_shipping_rate = rate; }
		void set_parent_county(County* in_county);
		void set_distancesNeighbors(std::multimap<double, Farm*>&);
		void set_neighborRadiusCalculated(const double);

		static void update_current_quarter_idx(int quarter_idx) { current_quarter_idx = quarter_idx; }
};

inline int Farm::get_id() const
{
	return id;
}
inline int Farm::get_cellID() const
{
	return cellID;
}
inline Farm_type* Farm::get_farm_type() const
{
  return farm_type;
}
inline double Farm::get_x() const
{
	return x_coordinate;
}
inline double Farm::get_y() const
{
	return y_coordinate;
}
inline double Farm::get_sus(unsigned int quarter_idx) const
{
	return sus.at(quarter_idx);
}
inline double Farm::get_inf(unsigned int quarter_idx) const
{
	return inf.at(quarter_idx);
}
inline double Farm::get_sus_max() const
{
    return sus_max;
}

inline double Farm::get_inf_max() const
{
    return inf_max;
}
inline std::string Farm::get_fips() const
{
	return fips;
}

inline County* Farm::get_parent_county() const
{
  return parent_county;
}
const inline std::multimap<double, Farm*>* Farm::get_distancesNeighbors()
{
	return &distancesNeighbors;
}
inline double Farm::get_neighborRadiusCalculated() const
{
	return neighborRadiusCalculated;
}


class Farm_type
{
public:
    Farm_type(std::string herd, std::vector<std::string> in_species, unsigned int index);
    ~Farm_type();
    std::string get_species() const { return species; }; //inlined
    unsigned int get_index() const { return index; }; //inlined
private:
    unsigned int index; // an index value relative to all other Farm_types
    std::string herd; // a binary code (ie '01') unique to this Farm_type, used to indicate species composition
    std::string species; //ie dairy, beef, name for this farm type

};

struct Prem_class
{
public:
    Prem_class(std::string tag, int idx) : tag(tag), idx(idx) {};
    std::string tag;
    int idx;
};

/// Derived class inheriting from Farm, containing additional info on replicate-specific statuses
/// and deleted at the end of each replicate. One object of this type is created for each
/// premises per simulation when it has any kind of status change (i.e. reporting, exposed,
/// prophylactic vaccination).
class Prem_status: public Farm
{
	private:
		std::string fileStatus; /// Current file status, one of: notDangerousContact, dangerousContact, reported
		std::string diseaseStatus; /// Current disease status, one of: exp, inf, imm for FMD simulations or btb for bTB simulations.
		std::vector<std::string> controlStatuses; ///< Vector of all control-related statuses ever effective for this premises
		std::string diagnosticStatus; /// Current diagnostic status, one of: started, complete, diagnostic type specific
		std::unordered_map<std::string, std::string> statusWhenWaitlisted; /// Record of premises status at the time of addition to control waitlist - the reason why control was targeted. Key is control type, value is status.
		std::unordered_map<std::string, std::string> statusWhenDiagnosticWaitlisted; /// Record of premises status at the time of addition to diagnostic waitlist - the reason why diagnostic was targeted. Key is diagnostic type, value is status.

		 std::unordered_map<std::string, int> testResult; ///< All diagnostic test results ever started for this premises, keyed by diagnostic type


		std::unordered_map<std::string, int> start; /// Start times for each status. Assumes unique names among file, disease, and control statuses due to mapping by status name.
		std::unordered_map<std::string, int> end; /// End times for each status. Assumes unique names among file, disease, and control statuses due to mapping by status name.

		std::vector<double> probPreventExposure; /// Probability of exposure being prevented, given effective control. One element for each control type currently effective.
		std::vector<double> probPreventTransmission; /// Probability of transmission being prevented, given effective control. One element for each control type currently effective.
		std::vector<Farm*> dangerousContacts; /// Farms that are dangerousContacts (set when premises is reported)
		std::unordered_map<Farm*, std::vector<bool>> potentialDCs; /// Map keyed by pointer to potentialDC, then by disease status, then bool if potentialDC is that status when this premises is reported

		// These vectors are "aligned" - each element corresponds to the same element in the others. Guessed this was easier and more flexible to search than tuple elements.
		std::vector<int> expTime;
		std::vector<Farm*> expSource;
		std::vector<int> expRoute;
		std::vector<std::string> expBlocked;
		std::vector<Shipment*> infectedShipmentsReceived; ///< Stores all the shipments with infected animals that this premises has received. Used for tracing. First element is oldest shipment.

		//For the FMD partial transition function. These are updated every time the FMD partial transition (PT) function is called
        double saved_evaluated_PTF_FMD = 0.0; //Everytime the PT function is evaluated it is saved here so it can be retrieved if we're still on the same timestep, instead of having to recalulate it multiple times each timestep with the same result.
        int saved_evaluated_PTF_FMD_timestep = -999; //This indicates which timestep the save d evaluated PT function is for. If this is the same as the current timestep, then we can reuse the saved value instead of recalculating.

        bool isVaccinated;
        std::unordered_map< std::string, double > currentUnvaccinatedPrevalence; ///< Proportion of animals of each type that are unvaccinated, keyed by types

        double get_inf_partial(int t, int quarter_idx, const Parameters* p,
                               const std::unordered_map<std::string, std::vector<double>>& normInf_map,
                               const std::unordered_map<std::string, int>& sp_counts);
        static double get_effective_n_infectious_animals(int time_since_infection,
                                                         int n_animals, const Parameters* p);

        //Btb-related
        std::vector<unsigned int> btb_infection_classes; //0=sus, 1=E1U, 2=E2U, 3=E2R, 4=IU, 5=IR.
        unsigned int btb_n_infectious; //IU + IR
        unsigned int btb_n_new_WH_infections = 0; ///< Tracks the number of new infections due to within-herd transmission each timestep.
        unsigned int btb_n_new_external_infections = 0; ///< Tracks the number of new infections from external sources (local spread) each timestep.
//        unsigned int btb_true_neg = 0;
//        unsigned int btb_false_pos = 0;
        void verify_premsize();

	public:
		Prem_status(Farm*);
		~Prem_status();

		void set_diagnosticStatus(std::string); //inlined

		void add_controlStatus(const std::string); //inlined
 		void add_exposureSource(int, Farm*, int, std::string); // confirmed record of source of infection, route, block
        void add_infected_shipment(Shipment* s) { infectedShipmentsReceived.push_back(s); };
 		void set_fileStatus(const std::string); //inlined
 		void set_diseaseStatus(const std::string); //inlined
 		void set_start(const std::string, const int); //inlined - set start time for status
 		void set_end(const std::string, const int); //inlined - set end time for status
 		void set_statusWhenWaitlisted(const std::string, const std::string);//inlined
 		void set_statusWhenDiagnosticWaitlisted(const std::string, const std::string);//inlined

		void add_testResult(const std::string, const int); // inlined
 		int get_testResult(std::string) const;
 		const std::vector<Shipment*>& get_received_inf_shipments() { return infectedShipmentsReceived; }

		double get_probPreventExposure() const;
 		void add_probPreventExposure(double); // inlined
 		void rem_probPreventExposure(double);

		double get_probPreventTransmission() const;
 		void add_probPreventTransmission(double); // inlined
		void rem_probPreventTransmission(double);

        bool get_isVaccinated() { return isVaccinated; }
        void vaccinate(double efficacy); //Determines the number of animals that are still susceptible after vaccination as  N_s~Bin(N_tot, 1-efficacy)
 		void unvaccinate(); //Removes effect of vaccination in the herd and sets N_s = N_tot.
 		std::unordered_map<std::string, int> get_currentSizeUnvaccinated();

 		std::string get_diseaseStatus() const; //inlined
 		std::string get_fileStatus() const; //inlined
 		int get_start(std::string) const;
 		int get_end(std::string) const;
 		std::vector<Farm*> get_dangerousContacts() const; //inlined
 		std::string get_diagnosticStatuses() const; //inlined
 		std::vector<std::string> get_controlStatuses() const; //inlined
 		bool is_exposureSource(Farm*);
 		void add_potentialDCInfo(Farm*, const std::vector<bool>&);//inlined
 		std::unordered_map<Farm*, std::vector<bool>>* get_potentialDCs();//inlined
 		void add_DC(Farm* f); //inlined
 		bool was_dcWhenWaitlisted(const std::string);

 		bool was_dcWhenDiagnosticWaitlisted(const std::string);

		bool is_onWaitlist(const std::string);
		bool is_onDiagnosticWaitlist(const std::string);

 		bool beenExposed() const; //inlined - not used?

        int when_infected(); //returns earliest exposure that was not blocked;

        //The following are only relevant for FMD runs with partial transition active.
        double get_inf_partial_as_unvaccinated(int t, int quarter_idx, const Parameters* p,
                                               const std::unordered_map<std::string, std::vector<double>>& normInf_map);
        double get_inf_partial_as_vaccinated(int t, int quarter_idx, const Parameters* p,
                                             const std::unordered_map<std::string, std::vector<double>>& normInf_map);
        //Btb-related
        int btb_seed();
//        void btb_initialize_false_pos();
        bool btb_update_size_change_quarter(); //Returns true if prem is still infected after update.
        double btb_get_prevalence_infectious();
        int btb_get_N_infected();
        unsigned int btb_get_new_WH_infections_since_last_t() { return btb_n_new_WH_infections; }
        unsigned int btb_get_new_external_infections_since_last_t() { return btb_n_new_external_infections; }
        void btb_expose_susceptibles(unsigned int n_to_expose);
        void btb_add_to_inf_classes(const std::vector<unsigned int>& n_in_btb_classes);
        void btb_subtract_from_inf_classes(size_t n, std::vector<unsigned int>& n_in_btb_classes);
        int btb_subtract_slaughter_shipment(size_t n_to_slaughter); //Subtracts n_to_slaughter animals from all infection classes based on prevalence and replace with susceptibles. Then return the number of the subtracted animals that are infected.
        void btb_remove_all_infected();
        bool btb_update_WH_spread(int quarter_idx); //Returns true if prem is infected.

        static double get_max_inf_partial(const Parameters* p,
                                          const std::unordered_map<std::string, double>& normInf_map,
                                          const std::unordered_map<std::string, int>& sp_counts);
};

inline std::string Prem_status::get_fileStatus() const
{
	return fileStatus;
}
inline std::string Prem_status::get_diseaseStatus() const
{
	return diseaseStatus;
}
inline std::vector<Farm*> Prem_status::get_dangerousContacts() const
{
	return dangerousContacts;
}
inline std::string Prem_status::get_diagnosticStatuses() const
{
	return diagnosticStatus;
}
inline std::vector<std::string> Prem_status::get_controlStatuses() const
{
	return controlStatuses;
}
inline void Prem_status::set_fileStatus(const std::string status)
{
	fileStatus = status;
}
inline void Prem_status::set_diseaseStatus(const std::string stat)
{
	diseaseStatus = stat;
}
inline void Prem_status::set_diagnosticStatus(std::string status)
{
	diagnosticStatus = status;
}
inline void Prem_status::add_controlStatus(const std::string status)
{
	controlStatuses.emplace_back(status);
}
inline void Prem_status::set_start(const std::string status, const int t)
{
	start[status] = t;
}
inline void Prem_status::set_end(const std::string status, const int t)
{
	end[status] = t;
}
inline void Prem_status::set_statusWhenWaitlisted(const std::string c_type, const std::string status)
{
	statusWhenWaitlisted[c_type] = status;
}

inline void Prem_status::set_statusWhenDiagnosticWaitlisted(const std::string d_type, const std::string status)
{
	statusWhenDiagnosticWaitlisted[d_type] = status;
}
inline void Prem_status::add_testResult(const std::string d_type, const int testres)
{
	testResult[d_type] = testres;
}
inline void Prem_status::add_probPreventExposure(double eff)
{
	probPreventExposure.emplace_back(eff);
}
inline void Prem_status::add_probPreventTransmission(double eff)
{
	probPreventTransmission.emplace_back(eff);
}
/// Checks whether or not a premises has been exposed to infection by checking for the presence of an "exp" start time
inline bool Prem_status::beenExposed() const
{
	return start.count("exp")==1;
}
inline std::unordered_map<Farm*, std::vector<bool>>* Prem_status::get_potentialDCs()
{
	return &potentialDCs;
}
inline void Prem_status::add_DC(Farm* f)
{
	dangerousContacts.emplace_back(f);
}
inline bool Prem_status::is_onWaitlist(const std::string c_type)
{
	return statusWhenWaitlisted.count(c_type)>0;
}
inline bool Prem_status::is_onDiagnosticWaitlist(const std::string d_type)
{
	return statusWhenDiagnosticWaitlisted.count(d_type)>0;
}

#endif //FARM_H
