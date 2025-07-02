//© 2024 Colorado State University
#ifndef Control_manager_h
#define Control_manager_h

#include <math.h>
#include "File_manager.h" // For parameters, Control_rule struct. Also includes shared functions.
#include "Grid_manager.h"
#include "County.h"
#include "State.h"
#include "Control_resource.h"

extern int verboseLevel;

class Farm;

/// Struct for a type of control (i.e. cull, vax), which can be implemented via one or more controlRules.
struct controlType
{
	std::string scale;
	std::tuple<double, double> effectiveness; /// Preventing exposure and preventing transmission
	std::string constraintType;
	std::tuple<double, double> implementDuration;
	std::tuple<double, double> effectiveDuration;
};

/// Struct of states of simulation variables to check for triggering control
struct statsForControlRules
{
	std::vector<Prem_status*>* newPremReports; ///< Vector of all farms reported in the last time step
	std::vector<Region_status*>* newCountyReports; ///< Vector of all counties reported in the last time step
	std::vector<Region_status*>* newStateReports; ///< Vector of all states reported in the last time step
	double percentIncrease; ///Percent increase in the number of farms reported from the previous timestep. 
	double decrease; ///Number of days without newPremReports or decreasing numbers of newPremReports.
	std::vector<std::string> controlsTriggered;
	std::unordered_map<Control_resource*, int> controlResourceLevels; // availability (keyed by control resource)
	std::unordered_map<std::string, double> waitlistFlag;
//	double waitlistFlag;
	int t; ///< Timestep
	// Could also include total number of premises, regions, animals, etc
};

/// Struct of premises/regions to be evaluated for waitlist removal, controlType of waitlist, and scale indicating premises or regions
struct waitlistGroup
{
	std::string controlType;
	std::string scaleType;
	std::vector<Farm*> premList;
	std::vector<std::string> regionList;
};

/// Stores controlTypes, controlRules, and initial states of Control_resources
class Control_manager
{
	private:
		int verbose; ///< Can be set to override global setting for console output
		const Parameters* p;
		Grid_manager* gridManager; // not const because neighbor calculations can change during simulations
		const std::unordered_map<int, Farm*>* allPrems;
		std::map<std::tuple<std::string, double, std::string>, bool> strategy; ///Map of all control triggers that turn off control. Keyed by control strategy. Value is whether control should be skipped or not (true/false).

		std::unordered_map<std::string, controlType*> allControlTypes; /// Key = name of control type, value = controlType struct
		std::unordered_map<std::string, std::unordered_map<std::string, Control_resource*>> controlResources; /// Map of all control resources - first key by control type, second key by ID (state, cellID, etc). Value is Control_resource
		std::unordered_map<std::string, std::unordered_map<int, int>> releaseSchedule; /// Map of all control resources boost timing - first key by control type, second key time third is amount of resource to add
		double daysIncreasing;
		void apply_rule(std::vector<Prem_status*>*, const controlRule&, std::vector<Farm*>&);
		void apply_regionRule(std::vector<Region_status*>*, const controlRule&, std::vector<std::string>&);
		void apply_percentRule(double percentIncrease, std::vector<Prem_status*>* reported, const controlRule& rule, std::vector<Farm*>& output);
		void apply_decreaseRule(const double& decrease, std::vector<std::string> controlsTriggered, const controlRule& rule);
		void apply_availabilityRule(std::vector<Prem_status*>* reported, std::unordered_map<Control_resource*, int> availability, std::vector<std::string> controlsTriggered, const controlRule& rule, std::unordered_map<std::string, double> flag, std::vector<Farm*>& output);
		void eval_controlRules(const controlRule& rule,const double& decrease, std::vector<std::string> controlsTriggered, std::vector<Prem_status*>* reported, std::unordered_map<Control_resource*, int> availability, std::unordered_map<std::string, double> flag,std::vector<Farm*>& output);
		void prioritize(std::string, std::vector<Farm*>&, std::vector<Farm*>&);
		void read_controlLocations(const Parameters*, std::string);
		void read_controlLocations_stateSum(const Parameters* p, std::string);

	public:
		Control_manager(const Parameters*, Grid_manager*);
		~Control_manager();

		

		const std::unordered_map<std::string, controlType*>* get_controlTypes() const;
		const std::unordered_map<std::string, std::unordered_map<std::string, Control_resource*>>* get_controlResources() const;
		const std::unordered_map<std::string, std::unordered_map<int, int>>* get_resourceBoostSchedule() const;
		std::vector<std::string> controlTriggered;
		void check_controlRules(statsForControlRules&, std::vector<waitlistGroup>&);
		void read_controlResourceBoost_nationalLimit(const Parameters* p, std::string controlType);
		void read_flexibleControlUpdateFile(const Parameters* p);		
		void filter_constraints(std::string, std::vector<Farm*>&, std::vector<Farm*>&,
			std::vector<Farm*>&, std::unordered_map<Control_resource*, int>&,
			std::unordered_map<std::string, std::unordered_map<Farm*, int>>&);
		void filter_constraints(std::string, std::vector<std::string>&,
			std::vector<std::string>&, std::vector<std::string>&,
			std::unordered_map<Control_resource*, int>*);

};

inline const std::unordered_map<std::string, controlType*>* Control_manager::get_controlTypes() const
{
	return &allControlTypes;
}

inline const std::unordered_map<std::string, std::unordered_map<std::string, Control_resource*>>* Control_manager::get_controlResources() const
{
	return &controlResources;
}

inline const std::unordered_map<std::string, std::unordered_map<int, int>>* Control_manager::get_resourceBoostSchedule() const
{
	return &releaseSchedule;
}

#endif // Control_manager_h
