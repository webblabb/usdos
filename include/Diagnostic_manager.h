#ifndef Diagnostic_manager_h
#define Diagnostic_manager_h

#include <math.h>
#include "File_manager.h" // For parameters, Diagnostic_rule struct. Also includes shared functions.
#include "Grid_manager.h"
#include "Farm.h"
#include "County.h"
#include "State.h"
#include "Diagnostic_resource.h"

extern int verboseLevel;

class Farm;

/// Struct for a type of Diagnostic (i.e. PCR, ELISA), which can be implemented via one or more diagnosticRules.
struct diagnosticType
{
	std::string scale;
	std::tuple<double, double> sensitivity; /// Accurately identifying infected premises
//	std::tuple<double, double> specificity; /// Accurately identifying non-infected premises
	std::string diagnosticConstraintType;
	std::tuple<double, double> startToCompleteDuration;
};

/// Struct of states of simulation variables to check for triggering diagnostics
struct statsForDiagnosticRules
{
	std::vector<Prem_status*>* newPremReports; ///< Vector of all farms reports (infected and reported) in the last time step
	std::vector<Prem_status*>* newPremSuspects; ///< Vector of all farms suspected (infected but not confirmed) in the last time step
	std::vector<Prem_status*>* newPremVaxs; ///< Vector of all farms that have effective vaccination in the last time step
	int t; ///< Timestep
	// Could also include total number of premises, regions, animals, etc
};

/// Struct of premises/regions to be evaluated for diagnostic waitlist removal, diagnosticType of waitlist, and scale indicating premises
struct diagnosticWaitlistGroup
{
	std::string diagnosticType;
	std::string scaleType;
	std::vector<Farm*> premList;
};

/// Stores diagnosticTypes, diagnosticRules, and initial states of Diagnostic_resources
class Diagnostic_manager
{
	private:
		int verbose; ///< Can be set to override global setting for console output
		const Parameters* p;
		Grid_manager* gridManager; // not const because neighbor calculations can change during simulations
		const std::unordered_map<int, Farm*>* allPrems;

		std::unordered_map<std::string, diagnosticType*> allDiagnosticTypes; /// Key = name of diagnostic type, value = diagnosticType struct
		std::unordered_map<std::string, std::unordered_map<std::string, Diagnostic_resource*>> diagnosticResources; /// Map of all diagnostic resources - first key by diagnostic type, diagnostic key by ID (state, cellID, etc). Value is Diagnostic_resource

		void apply_diagnosticRuleToSuspects(std::vector<Prem_status*>*, const diagnosticRule&, std::vector<Farm*>&);
		void apply_diagnosticRule(std::vector<Prem_status*>*, const diagnosticRule&, std::vector<Farm*>&);
		void apply_diagnosticRuleToVaxs(std::vector<Prem_status*>*, const diagnosticRule&, std::vector<Farm*>&);
		void prioritize(std::string, std::vector<Farm*>&, std::vector<Farm*>&);
	//	void read_diagnosticLocations(const Parameters*, std::string);

	public:
		Diagnostic_manager(const Parameters*, Grid_manager*);
		~Diagnostic_manager();

		const std::unordered_map<std::string, diagnosticType*>* get_diagnosticTypes() const;
		const std::unordered_map<std::string, std::unordered_map<std::string, Diagnostic_resource*>>* get_diagnosticResources() const;
		void check_diagnosticRules(statsForDiagnosticRules&, std::vector<diagnosticWaitlistGroup>&);
		void filter_diagnosticConstraints(std::string, std::vector<Farm*>&, std::vector<Farm*>&,
			std::vector<Farm*>&, std::unordered_map<Diagnostic_resource*, int>&,
			std::unordered_map<std::string, std::unordered_map<Farm*, int>>&);
		void filter_diagnosticConstraints(std::string, std::vector<std::string>&,
			std::vector<std::string>&, std::vector<std::string>&,
			std::unordered_map<Diagnostic_resource*, int>*);

};

inline const std::unordered_map<std::string, diagnosticType*>* Diagnostic_manager::get_diagnosticTypes() const
{
	return &allDiagnosticTypes;
}

inline const std::unordered_map<std::string, std::unordered_map<std::string, Diagnostic_resource*>>* Diagnostic_manager::get_diagnosticResources() const
{
	return &diagnosticResources;
}
#endif // Diagnostic_manager_h
