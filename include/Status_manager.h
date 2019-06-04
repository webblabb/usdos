#ifndef Status_manager_h
#define Status_manager_h

#include "shared_functions.h"
#include "Control_manager.h"
#include "Grid_manager.h"
#include "Shipment_manager.h"


#include <iterator> // for std::next
#include <utility> // for std::iter_swap

extern int verboseLevel;
class Farm;
class Region_status;
class Control_resource;

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
        Grid_manager* gridManager;

		const std::unordered_map<int, Farm*>* allPrems; ///< Pointer to Farms for interfacing with main
		const std::unordered_map<std::string, County*>* allCounties; ///< Pointer to Grid manager's county map
		const std::unordered_map<std::string, State*>* allStates; ///< Pointer to Grid manager's state map
		const std::unordered_map<std::string, controlType*>* allControlTypes; ///< Pointer to Control manager's control types
		const std::unordered_map<std::string, std::unordered_map<std::string, Control_resource*>>* controlResources; ///< Pointer to Control manager's control types
		const std::unordered_map<std::string, std::unordered_map<int, int>>* controlReleaseSchedule; ///< Pointer to Control manager's control types

		unsigned int recentNotSus; ///< Placeholder for last farm that became not-susceptible during the last timestep
		std::tuple<double, double> pastEndTime; ///< A time lag used for indefinite or permanent statuses
		int nPrems; ///< Total number of premises, needed to return number of susceptibles

		int verbose; ///< Can be set to override global setting for console output

		int firstRepTime; //to record first reported time

		std::unordered_map<std::string, statusShift> statusSequences; ///< Statuses that auto-advance to next status and duration of next status, includes disease, file, and control statuses

		std::unordered_map<std::string, statusList<Prem_status*>> fileStatuses; ///< Vectors of farms with particular file statuses, with current validity placeholder
		std::unordered_map<std::string, statusList<Prem_status*>> diseaseStatuses; ///< Vectors of farms with particular disease statuses, with current validity placeholder
		std::unordered_map<std::string, statusList<Prem_status*>> controlStatuses; ///< Vectors of farms with particular control statuses, with current validity placeholder
		std::unordered_map<std::string, statusList<Region_status*>> regionControlStatuses; ///< Vectors of regions with particular control statuses, with current validity placeholder - region statuses only need to be updated for control, disease and file statuses tracked at premises level

		std::vector<std::string> species;

		std::unordered_map<int, Prem_status*> changedStatus; ///< Premises that have had any status change
		std::unordered_map<std::string, Region_status*> changedCoStatus; ///< Counties that have had any status change
		std::unordered_map<std::string, Region_status*> changedStateStatus; ///< States that have had any status change
		std::vector<Prem_status*> newPremReports; ///< Starting index for newly reported premises
		std::vector<Region_status*> newCoReports; ///< Starting index for newly reported counties
		std::vector<Region_status*> newStateReports; ///< Starting index for newly reported states

		std::vector<Farm*> notSus; ///< Farms that are in any disease state except susceptible (are not eligible for local spread exposure)

		std::vector<std::tuple<Farm*, Farm*, int, std::string>> sources; ///< Exposed farm, source of infection, type of spread (0=local, 1=ship), controlPrevented string ("shipBan" or "premControl")
 		std::vector<std::string> reportedCounties;
 		std::vector<std::string> reportedStates;
		std::vector<std::tuple<Farm*, Farm*, int, double>> exposureForEval; ///< Exposures to be confirmed against control in this timestep(destination, origin, route, probability of exposure)
		std::unordered_map<std::string, std::vector<Farm*>> waitlist; ///< Waitlists of premises for each control type
		std::unordered_map<std::string, std::vector<std::string>> waitlistRegion; ///< Waitlists of regions for each control type
		std::vector<int> dcsPerIP;

		std::unordered_map<Control_resource*, int> controlResourceLevels; // availability (keyed by control resource)
		std::unordered_map<std::string, std::unordered_map<Farm*, int>> partiallyControlledPrems; // key control type, then farm, value is animals remaining to be controlled

		int verify_premStatus(Farm*);
		void get_seedCos(std::vector<std::string>&);
		void set_status(Farm*, int, std::string, std::tuple<double,double>, std::tuple<double, double> controlEffect = std::make_tuple(0,0));
		void set_regionStatus(std::string id, std::string regionType, int, std::string, std::tuple<double,double>, std::tuple<double, double> controlEffect = std::make_tuple(0,0));
		void report_countyAndState(Farm* f, int t);
		void eval_premExpPrevByVaccination(Farm* ofarm, Farm* dfarm, double trueP, int t, bool& transPrevented, bool& expPrevented);
		bool eval_premTransmission(Farm*);
		bool eval_premExposure(Farm*);
		void expose(std::vector<std::pair<Farm*, int>>&, int); //Takes a vector of pairs, each pair is the farm to be exposed and the specific latency to be used for this particular exposure.
		void update(int t, std::string, statusList<Prem_status*>&);
		void update(int t, std::string, statusList<Region_status*>&);
		void add_premSource(int, Farm*, Farm*, int, std::string); //inlined
		int count_allDCPrems(const std::string);

	public:
		Status_manager(std::vector<Farm*>&, const Parameters*, Grid_manager*,
		  Control_manager*);
		~Status_manager();

		void add_premForEval(Farm*, Farm*, int, double); //inlined
		void updateDisease(int t);
		void updateControl(int t);

		void eval_exposure(int); // check for control before exposure
		void filter_shipments(std::vector<Shipment*>&, int); // check shipBans, recipient disease statuses
		void get_premsWithStatus(std::vector<std::string> status_vector, std::vector<Farm*>& output); // returns vector of Prem_status*s with the statuses provided in the argument status vector
		void get_premsWithStatus(std::string, std::vector<Farm*>&); // returns vector of Prem_status*s with status
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
		std::string formatRepSummary(int, int, double);
		std::string formatDetails(int, int);

		void add_waitlistMembers(int);
		void update_ControlResources(int t);
		void add_implemented(int);
		void add_potentialDC(Farm*, Farm*, std::unordered_map<std::string, bool>&);

        Prem_status* get_correspondingPremStatus(int fid);
        const std::unordered_map<std::string, double>& get_normInf_map() const;
        const std::unordered_map<std::string, double>& get_normSus_map() const;
};

/// Temporary list of exposure sources needing evaluation of premises-level control. Used by
/// Grid_checker::stepThroughCells and Status_manager::filter_shipments
inline void Status_manager::add_premForEval(Farm* toBeExposed, Farm* exposedBy, int route, double trueP)
{
	exposureForEval.emplace_back(std::make_tuple(toBeExposed, exposedBy, route, trueP));
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
