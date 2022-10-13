#ifndef Status_manager_h
#define Status_manager_h

#include "shared_functions.h"
#include "Control_manager.h"
#include "Grid_manager.h"
#include "Shipment_manager.h"
#include "Diagnostic_manager.h"


#include <iterator> // for std::next
#include <utility> // for std::iter_swap

extern int verboseLevel;
class Farm;
class Region_status;
class Control_resource;
class Diagnostic_resource;

/// Struct containing vector of all farms/regions that have/had a given status, with placeholders to indicate current validity
template<typename T>
struct statusList
{
	std::vector<T> units;///< List of all farms/regions with this status, in chronological order of start time
	unsigned int lo; ///< Before this placeholder, end times are less than t (have already passed)
	unsigned int hi; ///< After this placeholder, start times are greater than t (have not started yet)
};

/// Struct defining a transition from one status to the next. A statusShift exists for a status
/// A that transitions to status B, where 'duration' is the mean and variance of time in days for
/// A to transition to B, and 'next' is status B. When status A expires, Status_manager::updates()
/// assigns the status in 'next'
struct statusShift
{
	std::tuple<double,double> duration; ///< Mean and variance of duration of status
	std::string next; ///< Following status
};

/// 	Keeps track of disease and control statuses during a simulation.
/// 	Filters exposures from local and shipment spread with any control measures in place.
class Status_manager
{
	private:
		std::vector<Farm*> seededFarms;
		const Parameters* parameters;
		Control_manager* controlManager;
		Diagnostic_manager* diagnosticManager;
        Grid_manager* gridManager;

		const std::unordered_map<int, Farm*>* allPrems; ///< Pointer to Farms for interfacing with main
		const std::unordered_map<std::string, County*>* allCounties; ///< Pointer to Grid manager's county map
		const std::unordered_map<std::string, State*>* allStates; ///< Pointer to Grid manager's state map
		const std::unordered_map<std::string, controlType*>* allControlTypes; ///< Pointer to Control manager's control types
		const std::unordered_map<std::string, std::unordered_map<std::string, Control_resource*>>* controlResources; ///< Pointer to Control manager's control types
		const std::unordered_map<std::string, std::unordered_map<int, int>>* controlReleaseSchedule; ///< Pointer to Control manager's control types
		const std::unordered_map<std::string, diagnosticType*>* allDiagnosticTypes; ///< Pointer to Diagnostic manager's diagnostic types
		const std::unordered_map<std::string, std::unordered_map<std::string, Diagnostic_resource*>>* diagnosticResources; ///< Pointer to Diagnostic manager's control types

		unsigned int recentNotSus; ///< Placeholder for last farm that became not-susceptible during the last timestep
		std::tuple<double, double> pastEndTime; ///< A time lag used for indefinite or permanent statuses
		int nPrems; ///< Total number of premises, needed to return number of susceptibles

		int verbose; ///< Can be set to override global setting for console output

		int firstRepTime; //to record first reported time

		std::unordered_map<std::string, statusShift> statusSequences; ///< Statuses that auto-advance to next status and duration of next status, includes disease, file, and control statuses

		std::unordered_map<std::string, statusList<Prem_status*>> fileStatuses; ///< Vectors of farms with particular file statuses, with current validity placeholder
		std::unordered_map<std::string, statusList<Prem_status*>> diseaseStatuses; ///< Vectors of farms with particular disease statuses, with current validity placeholder
		std::unordered_map<std::string, statusList<Prem_status*>> diagnosticStatuses; ///< Vectors of farms with particular diagnostic statuses, with current validity placeholder
	//	std::unordered_map<std::string, statusList<Region_status*>> regionDiagnosticStatuses; ///< Vectors of regions with particular diagnostic statuses, with current validity placeholder - region statuses only need to be updated for diagnostic, disease and file statuses tracked at premises level
		std::unordered_map<std::string, statusList<Prem_status*>> controlStatuses; ///< Vectors of farms with particular control statuses, with current validity placeholder
		std::unordered_map<std::string, statusList<Region_status*>> regionControlStatuses; ///< Vectors of regions with particular control statuses, with current validity placeholder - region statuses only need to be updated for control, disease and file statuses tracked at premises level

		std::vector<std::string> species;

		std::unordered_map<int, Prem_status*> changedStatus; ///< Premises that have had any status change
		std::unordered_map<std::string, Region_status*> changedCoStatus; ///< Counties that have had any status change
		std::unordered_map<std::string, Region_status*> changedStateStatus; ///< States that have had any status change
		std::vector<Prem_status*> newPremReports; ///< Starting index for newly reported premises
		std::vector<Region_status*> newCoReports; ///< Starting index for newly reported counties
		std::vector<Region_status*> newStateReports; ///< Starting index for newly reported states

		std::vector<Prem_status*> newPremSuspects; ///< Starting index for newly suspects premises
		std::vector<Prem_status*> newPremVaxs; ///< Starting index for newly effective vaccinated premises

		std::vector<Farm*> notSus; ///< Farms that are in any disease state except susceptible (are not eligible for local spread exposure)

		std::vector<std::tuple<Farm*, Farm*, int, std::string>> sources; ///< Exposed farm, source of infection, type of spread (0=local, 1=ship), controlPrevented string ("shipBan" or "premControl")
 		std::vector<std::string> reportedCounties;
 		std::vector<std::string> reportedStates;
		std::vector<std::tuple<Farm*, Farm*, int, double, Shipment*>> exposureForEval; ///< Exposures to be confirmed against control in this timestep(destination, origin, route, probability of exposure)
		std::unordered_map<std::string, std::vector<Farm*>> waitlist; ///< Waitlists of premises for each control type
		std::unordered_map<std::string, std::vector<std::string>> waitlistRegion; ///< Waitlists of regions for each control type
		std::unordered_map<std::string, std::vector<Farm*>> diagnosticWaitlist; ///< Waitlists of premises for each diagnostic type
		std::unordered_map<std::string, std::vector<Farm*>> toTest; ///< premises to test for each diagnostic type
		std::vector<int> dcsPerIP;

		std::unordered_map<Control_resource*, int> controlResourceLevels; // availability (keyed by control resource)
		std::unordered_map<std::string, std::unordered_map<Farm*, int>> partiallyControlledPrems; // key control type, then farm, value is animals remaining to be controlled
		std::unordered_map<Diagnostic_resource*, int> diagnosticResourceLevels; // availability (keyed by diagnostic resource)
		std::unordered_map<std::string, std::unordered_map<Farm*, int>> partiallyTestedPrems; // key diagnostic type, then farm, value is animals remaining to be tested
	    std::vector<std::string> speciesOnPrems; ///< List of species on all farms provided in premises file

	    std::string btb_slaughter_ofn; ///< Filename for storing slaughter shipment information.
	    std::ofstream btb_slaughter_ofs; ///< File stream for storing slaughter shipment information.
        std::unordered_set<Prem_status*> btb_all_infected; ///< Keeps track of the number of unique premises that have had infection.
        std::unordered_set<Prem_status*> btb_all_recovered; ///< Keeps track of the number of unique premises that have recovered.
        std::unordered_set<Prem_status*> btb_all_reported; ///< Keeps track of the number of unique premises that have been reported.
        int btb_n_freshly_infected_animals = 0; ///< Keeps track of the number of unique infection events of animals.

		int verify_premStatus(Farm*);
		void get_seedCos(std::vector<std::string>&);
		void set_status(Farm*, int, std::string, std::tuple<double,double>, std::tuple<double, double> controlEffect = std::make_tuple(0,0));
		void set_regionStatus(std::string id, std::string regionType, int, std::string, std::tuple<double,double>, std::tuple<double, double> controlEffect = std::make_tuple(0,0));
		void report_countyAndState(Farm* f, int t);
		void eval_premExpPrevByVaccination(Farm* ofarm, Farm* dfarm, double trueP, int t, int quarter_idx, bool& transPrevented, bool& expPrevented);
		bool eval_premTransmission(Farm*);
		bool eval_premExposure(Farm*);
		void expose(std::vector<std::pair<Farm*, int>>&, int); //Takes a vector of pairs, each pair is the farm to be exposed and the specific latency to be used for this particular exposure.
		void update(int t, std::string, statusList<Prem_status*>&);
		void update(int t, std::string, statusList<Region_status*>&);
		void add_premSource(int, Farm*, Farm*, int, std::string); //inlined
		int count_allDCPrems(const std::string);
		int count_allDDCPrems(const std::string);

		void btb_updateWithinHerd(int t, int quarter_idx);
		Prem_status* btb_make_trace(Prem_status* pst, int current_depth = 0);

	public:
		Status_manager(std::vector<Farm*>&, const Parameters*, Grid_manager*,
		  Control_manager*, Diagnostic_manager*);
		~Status_manager();

		void add_premForEval(Farm*, Farm*, int, double, Shipment* s=nullptr); //inlined
		void updateDisease(int t);
		void updateFile(int t);
		void updateControl(int t);
		void updateDiagnostics(int t);

		bool eval_probHerdPositive(Farm* f, int t, int quarter_idx, std::string); //evaluate diagnostic test results
		void eval_testResult(int t, int quarter_idx); // check for diagnostics before determining next status
		void eval_slaughterSurveillance(int t, int quarter_idx, int replicate, std::vector<Farm*> prems);
		void wipe_markets(); //Removes all infected animals from markets and replaces them with susceptibles. Only ever occurs if set by config option 23.

		void eval_exposure(int t, int quarter_idx); // check for control before exposure
		void filter_shipments(std::vector<Shipment*>&, int); // check shipBans, recipient disease statuses
		void get_premsWithStatus(std::vector<std::string> status_vector, std::vector<Farm*>& output); // returns vector of Farm*s with the statuses provided in the argument status vector
		void get_premsWithStatus(std::string, std::vector<Farm*>&); // returns vector of Farm*s with status
		void get_premStatusesWithStatus(std::vector<std::string> status_vector, std::vector<Prem_status*>& output); //Returns vector of Prem_status*s with the statuses provided in the argument status vector
		void get_premStatusesWithStatus(std::string status, std::vector<Prem_status*>& output); //Returns vector of Prem_status*s with the statuses provided by the argument status
		int get_totalPremsWithStatus(std::string); // get number of premises _ever_ with this status
		int numPremsWithStatus(std::string); // get number of Prem_status*s with disease status
		int numPremsWithFileStatus(std::string s); // get number of Prem_status*s with file status
		int get_totalPremsWithFileStatus(std::string); // get number of premises _ever_ with this file status status
		int get_numRegionsWithControlStatus(std::string, std::string);
 		int get_numCountiesReported() const; //inlined
 		int get_numStatesReported() const; //inlined
		void newNotSus(std::vector<Farm*>&); //inlined
		std::string getAny_fileStatus(Farm*) const;
		std::string getAny_diseaseStatus(Farm*) const; // exists to output "sus" in case of no Prem_status
		std::string getAny_diagnosticStatus(Farm*);
		std::string formatRepSummary(int, int, double);
		std::string formatDetails(int, int);

		void add_waitlistMembers(int);
		void add_diagnosticWaitlistMembers(int);
		void update_ControlResources(int t);
		void add_implemented(int);
		void add_started(int);
		void add_potentialDC(Farm*, Farm*, std::vector<bool>&);

        Prem_status* get_correspondingPremStatus(int fid);
        const std::unordered_map<std::string, std::vector<double>>& get_normInf_map() const;
        const std::unordered_map<std::string, std::vector<double>>& get_normSus_map() const;

        //Members specific to btb only.
        void btb_updateDisease(int t, int quarter_idx);
        void btb_updatePremSizeChangesQuarter(int t); //Goes through all Prem_status objects that have exposed animals and updates the prevalence according to the current size of the prem. Used after a change of quarter gives all premises new sizes.
        void btb_evalExposure(int t); //Takes a vector of shipments and goes through them, sampling animals from the various btb infection classes and adds them to the shipment and then updates the prevalence for the sending and receiving premises.
        void btb_updateRecords(); ///< Updates the records of how many unique premises and farms have been affected.
};

/// Temporary list of exposure sources needing evaluation of premises-level control. Used by
/// Grid_checker::stepThroughCells and Status_manager::filter_shipments
inline void Status_manager::add_premForEval(Farm* toBeExposed, Farm* exposedBy, int route,
                                            double trueP, Shipment* s)
{
	exposureForEval.emplace_back(std::make_tuple(toBeExposed, exposedBy, route, trueP, s));
}
inline int Status_manager::get_numCountiesReported() const
{
	return reportedCounties.size();
}
inline int Status_manager::get_numStatesReported() const
{
	return reportedStates.size();
}

#endif // Status_manager_h
