#ifndef file_manager_h
#define file_manager_h

#include <fstream>
#include "shared_functions.h"
#include "Local_spread.h"
// shared_functions includes iostream, sstream, string, vector

extern int verboseLevel;

/// Struct for a rule relating to how controls are implemented. If (trigger>threshold),
/// apply action to target, in priority order. For distance-based control, a priority
/// value of "closest" may result in some additional computation time, as any entire grid
/// cell entirely within-radius has all farms designated as in-radius without distance
/// calculations. However, prioritizing by distance means those distances need to actually
/// be calculated.
struct controlRule
{
	std::string trigger;
	double threshold;
	std::string action;
	double target;
	std::string priority;
	double radiusSquared; ///< Used only if target = radius of neighbors, stored to avoid multiple calculations
};

/// Struct containing all parameters loaded from file
struct Parameters
{
	// output parameters
	std::string batch; ///< Prefix for summary and detailed output files
	int printSummary;
	int printDetail;
	int printCells;
	int printShipments;
	int printControl;

	// general parameters
	std::string premFile; ///< File containing tab-delimited premises data: ID, FIPS, x, y, population
	std::string fipsFile;
	std::vector<std::string> species;
	int timesteps;
	bool useMaxPrems;
	int maxInfectiousPrems;
	int start_day_option;
	int start_day;
	int replicates;
	int verboseLevel;
	bool pairwiseOn;
	bool reverseXY;

	// infection parameters
	std::string seedSource;
	std::string dataKernelFile; ///< Name of file containing local spread probabilities by distance (set kernelType to 1)
	std::string seedSourceType;
	std::unordered_map<std::string,double> susExponents; ///< Species-specific exponents for susceptibility (q in USDOSv1)
	std::unordered_map<std::string,double> infExponents; ///< Species-specific exponents for infectiousness (p in USDOSv1)
	std::unordered_map<std::string,double> susConsts; ///< Species-specific constants for susceptibility (S in Tildesley ProcB 2008)
	std::unordered_map<std::string,double> infConsts; ///< Species-specific constants for infectiousness (T in Tildesley ProcB 2008)
	std::tuple<double,double> latencyParams;
	std::tuple<double,double> infectiousParams;


    //partial transistion flag
    int partial;
    // partial transition parameters
    std::vector<double> partialParams;
     // Species-specific birth rates
	//std::unordered_map<std::string,double> popParams;

	// local spread kernel object (for grid manager and checker)
	int kernelType;
	std::vector<double> kernelParams;
	Local_spread* kernel;

	// grid parameters
	std::string cellFile;
	std::vector<int> densityParams;
	int uniformSide;

	// shipment parameters
	bool shipments_on;
	std::string shipment_kernel;
	int usamm_version;
//	std::vector<int> shipMethods;
    int shipMethods;
	std::vector<int> shipMethodTimeStarts;
	int shipPremAssignment;
	std::vector<std::string> USAMM_parameter_files;
	std::vector<std::string> USAMM_temporal_order;
	std::vector<int> USAMM_temporal_start_points;
	std::vector<int> USAMM_temporal_n_days;
	std::vector<std::string> USAMM_ocov_files;
	std::vector<std::string> USAMM_dcov_files;
	std::vector<std::string> USAMM_supernode_files;
	bool exposed_shipments;
	std::vector<std::string> statuses_to_generate_shipments_from; //These are the statuses to be considered when generating shipments.

	// control parameters - maps keyed by controlType name
	bool control_on;
	std::vector<std::string> controlTypes;
	std::unordered_map<std::string, std::string> controlScales;
	std::unordered_map<std::string, std::string> constraintFunctions;
	std::unordered_map<std::string, std::vector<double>> constraintFuncParams;
	std::unordered_map<std::string, std::vector<std::string>> constraintFuncFiles;
	std::unordered_map<std::string, std::vector<std::string>> constraintFuncFileTypes;
	std::unordered_map<std::string, std::tuple<double,double>> implementToEffectiveLag;
	std::unordered_map<std::string, std::tuple<double,double>> effectiveToInactiveLag;
	std::unordered_map<std::string, std::tuple<double,double>> effectiveness;

	// control rules
	std::vector<controlRule> controlRules;
	std::vector<std::string> dcControlTypes;

	// reporting parameters
	std::tuple<double,double> indexReportLag;
	std::tuple<double,double> nonDCReportLag;
	std::tuple<double,double> dcReportLag;

	// Dangerous Contacts parameters
	bool dangerousContacts_on;
	std::unordered_map<std::string, double> dcRiskScale;
	double maxDCScale;
};

/// Loads and checks parameters from configuration file
class File_manager
{
	private:
		int verbose; ///< Can be set to override global setting for console output
		std::vector<std::string> pv; ///< Parameter vector for reading in from file
		Parameters params;

		bool checkMeanVar(std::string&, int, std::string);
		bool checkPositive(std::vector<int>&, int);
		bool checkPositive(std::vector<double>&, int);
		bool checkZeroToOne(std::vector<double>&, int);
		bool limitedValues(std::vector<std::string>&, std::vector<std::string>&, int);

	public:
		File_manager();
		~File_manager();
		const Parameters* getParams(); // inlined
		const std::string getSettings(std::string&);
		void readConfig(std::string&);
};

inline const Parameters* File_manager::getParams()
{
	return &params;
}

#endif // File_manager_h
