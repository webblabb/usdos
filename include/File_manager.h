//© 2024 Colorado State University
#ifndef file_manager_h
#define file_manager_h

#include <fstream>
#include "shared_functions.h"
#include "Local_spread.h"
// shared_functions includes iostream, sstream, string, vector

extern int verboseLevel;
enum class InfectionType { FMD = 0, BTB };
enum class DiseaseStatus { SUS = 0, EXP, INF, IMM, BTB, REC };

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
	double delay; ///< This rule is only used during flexible control strategy runs and is initialized with a value only if a flex.txt file is provided.
};

struct BtbWithinHerdParams
{
    //Number of animals to start with in infection classes.
    int S_0 = 750; //Susceptible
    int E1U_0 = 0; //Exposed, unreactive, no immune response
    int E2R_0 = 0; //Exposed, reactive, immune response
    int E2U_0 = 0; //Exposed, unreactive, immune response
    int IU_0 = 0; //Infectious, unreactive
    int IR_0 = 0; //Infectious, reactive

    //Premises-level parameters
    double export_rate = 0.9; // mu
    double export_prop_q4 = 0.4; // q4_mu or q_mu in paper, proportion of export rate that is contributed by 4th quarter
    double export_rate_q1 = export_rate * (1 - export_prop_q4) / 3.0; // mu_q1
    double export_rate_q2 = export_rate * (1 - export_prop_q4) / 3.0; // mu_q2
    double export_rate_q3 = export_rate * (1 - export_prop_q4) / 3.0; // mu_q3
    double export_rate_q4 = export_rate * export_prop_q4; // mu_q4
    std::vector<double> export_rate_v = {
        export_rate_q1,
        export_rate_q2,
        export_rate_q3,
        export_rate_q4
    };
    double dairy_birth_rate = 0.05;
    double dairy_mortality_rate = 0.05; // nu or eta in paper
    double beef_mortality_rate = 0.02; // nu or eta in paper
    double beef_import_rate = 0.49; //alpha
    double dairy_import_rate = 0.09; //alpha

    //Transmision parameters
    double transmission_rate = 25.0; // beta?
    double cattle_cattle_contact_rate = 0.3; // p1
    double cattle_wildl_contact_rate_q1 = 0.0; // p2_q1
    double cattle_wildl_contact_rate_q2 = 0.0; // p2_q2
    double cattle_wildl_contact_rate_q3 = 0.0; // p2_q3
    double cattle_wildl_contact_rate_q4 = 0.0; // p2_q4
    std::vector<double> contact_rate_v = {
        cattle_wildl_contact_rate_q1,
        cattle_wildl_contact_rate_q2,
        cattle_wildl_contact_rate_q3,
        cattle_wildl_contact_rate_q4
    };
//    double cattle_fomite_contact_rate = 0.1; // p3
    double cattle_fomite_contact_rate = 0.0; // p3
    double prop_infected_imports = 0.0; // pi
    double prop_reactive_are_infectious = 0.5; // psi_1
    double prop_unreactive_are_infectious = 0.5; // psi_2
    double prop_unreactive_exposed_wo_immune_resp = 0.5; // psi_3

    //Latency parameters



    //sigma_1 is the rate of transition from class E1 to E2 or E3 (mounting an immune response) and is drawn from a gamma distribution
    //with these parameters
    double sigma_1_mean = 14.0;// sigma1
    double sigma_1_rate = 0.5; // sigma1_rate
    double sigma_1_scale = 1.0 / sigma_1_rate; //Just for use with scipys gamma rv fun that takes scale rather than rate.
    double sigma_1_shape = sigma_1_mean * sigma_1_rate; // sigma1_rate

    //sigma_2 is the rate of transition from exposed and reactive (E2) to infectious (I1, I2) and is drawn from a gamma distribution
    //with these parameters
    double sigma_2_mean = 12.0; // sigma2
    double sigma_2_rate = 1.7; // sigma2_rate
    double sigma_2_scale = 1.0 / sigma_2_rate;
    double sigma_2_shape = sigma_2_mean * sigma_2_rate; // sigma2_shape

    //sigma_2_star is the rate of transition from exposed and unreactive (E3) to infectious and unreactive (I1) and is lower than sigma_2 (takes longer for unreactive to become infectious). iIt  is drawn from a gamma distribution
    double sigma_2_star_rate = sigma_2_rate; // sigma2_star_rate
    double sigma_2_star_scale = 1.0 / sigma_2_star_rate;

    double phi_sigma_2 = 1.6; // phi
    double sigma_2_star_shape = sigma_2_shape * phi_sigma_2; // sigma2_star_shape

    //delta_1 is the rate at which an unreactive exposed animal become reactive
    double delta_1_rate = sigma_1_rate; // delta1_rate
    double delta_1_scale = 1.0 / delta_1_rate;

    double phi_delta_1 = 1.6; // phi_delta1
    double delta_1_shape = sigma_1_shape * phi_delta_1; // delta1_shape

    //delta_2 is the rate at which an infectious unreactive animal becomes reactive
    double delta_2_rate = sigma_2_rate; // delta2_rate
    double delta_2_scale = 1.0 / delta_2_rate;

    double phi_delta_2 = 0.7;// phi_delta2
    double delta_2_shape = sigma_2_shape*phi_delta_2; // delta2_shape


    //Diagnostic testing parameters
    double response_prob_1 = 0.0; // pa_1, rho_1 Probability of an exposed animal without immune response (E1U) that mounts an immune response (leaves E1U) becomes reactive to testing (goes into E2R).
    double response_prob_2 = 1.0; // pa_2, rho_2 Probability of an unreactive exposed animal with immune response (E2U) becomes reactive to testing, meaning it moves into E2R.
    double response_prob_3 = 0.0; // pa_3, rho_3 Probability that an infectious animal that is not reactive to testing (IU) becomes reactive (goes into IR).
    double test_sensitivity = 0.9; // t1
    double test_specificity = 0.0; // t2
    double specificity_scaling_factor = 0.0; // pa_t2, for sens. analysis?

    std::map<int, int> month_to_quarter = { {1, 0}, {2, 0}, {3, 0}, {4, 1}, {5, 1}, {6, 1},
                                            {7, 2}, {8, 2}, {9, 2}, {10, 3}, {11, 3}, {12, 3} };

    //Must be called after finishing setting of parameters in File_manager.
    void update_latency_parameters()
    {
        sigma_1_scale = 1.0 / sigma_1_rate;
        sigma_1_shape = sigma_1_mean * sigma_1_rate;
        sigma_2_scale = 1.0 / sigma_2_rate;
        sigma_2_shape = sigma_2_mean * sigma_2_rate;
        sigma_2_star_rate = sigma_2_rate;
        sigma_2_star_scale = 1.0 / sigma_2_star_rate;
        sigma_2_star_shape = sigma_2_shape * phi_sigma_2;
        delta_1_rate = sigma_1_rate;
        delta_1_scale = 1.0 / delta_1_rate;
        delta_1_shape = sigma_1_shape * phi_delta_1;
        delta_2_rate = sigma_2_rate;
        delta_2_scale = 1.0 / delta_2_rate;
        delta_2_shape = sigma_2_shape*phi_delta_2;
    }


};

/// Struct for a rule relating to how diagnostics are implemented. If (trigger>threshold),
/// apply action to target, in priority order. For distance-based control, a priority
/// value of "closest" may result in some additional computation time, as any entire grid
/// cell entirely within-radius has all farms designated as in-radius without distance
/// calculations. However, prioritizing by distance means those distances need to actually
/// be calculated.
struct diagnosticRule
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
    //"Global" constants.
    const int NUM_INF_TYPES = 2; //FMD, BTB
    const static int NUM_DISEASE_STATUSES = 6; //SUS, EXP, INF, IMM, BTB, REC
    const static std::map<DiseaseStatus, std::string> DiseaseStatusToString;
    const static std::map<std::string, DiseaseStatus> StringToDiseaseStatus;
    const std::vector<int> DAYS_PER_MONTH = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

	// output parameters
	std::string batch; ///< Prefix for summary and detailed output files
	std::string batchDateTime;
	int printSummary;
	int printDetail;
	int printCells;
	int printShipments;
	int printControl;
	int printDiagnostic;
	bool generate_shipment_network;

	// general parameters
	std::string premFile; ///< File containing tab-delimited premises data: ID, FIPS, x, y, population
	std::string fipsFile;
	std::vector<std::string> species;
	int timesteps;
	int days_per_timestep = 1;
	bool useMaxPrems;
	int maxInfectiousPrems;
	int start_day_option;
	int start_day;
	int replicates;
	int verboseLevel;
	bool pairwiseOn;
	bool reverseXY;
	InfectionType infectionType = InfectionType::FMD;

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
    BtbWithinHerdParams btbWithinHerdParams;
    bool within_spread_at_markets; //Controls whether within-premises spread takes place at markets.
    double market_shipment_exposure_prob = 1.0;
    bool wipe_markets_each_t = false; //Controls wheter markets should be emptied of infected animals at the end of each time step.
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
    double shipment_rate_factor = 1.0; //Scales up the amount of shipments linearly by this factor. Applied directly on USAMM parameters N, lambda or c depending on USAMM version when fetching parameter from the USAMMvX_parameters object.
	std::vector<int> shipMethodTimeStarts;
	std::vector<std::string> USAMM_parameter_files;
	std::vector<unsigned int> USAMM_post_lengths; //Number of lines in the posterior files, saved here to improve speed by not having to check each time a new sample is made.
	std::vector<std::string> USAMM_temporal_order;
	std::vector<int> USAMM_temporal_start_points;
	std::vector<int> USAMM_temporal_n_days;
	std::vector<int> USAMM_temporal_n_timesteps;
	std::vector<std::string> USAMM_ocov_files;
	std::vector<std::string> USAMM_dcov_files;
	std::vector<std::string> USAMM_supernode_files;
	std::vector<std::string> statuses_to_generate_shipments_from; //These are the statuses to be considered when generating shipments.
	std::vector<std::string> statuses_to_generate_slaughter_shipments_from; //These are the statuses to be considered when generating slaughter shipments.

	// control parameters - maps keyed by controlType name
	bool control_on;
	std::vector<std::string> controlTypes;
	std::string flexibleControl;
	std::string flexibleControlFiles;
	std::unordered_map<std::string, std::string> controlScales;
	std::unordered_map<std::string, std::string> constraintFunctions;
	std::unordered_map<std::string, std::vector<double>> constraintFuncParams;
	std::unordered_map<std::string, std::vector<std::string>> constraintFuncFiles;
	//std::unordered_map<std::string, std::vector<std::string>> flexibleControlFiles; //flexible control files
	std::unordered_map<std::string, std::vector<std::string>> constraintFuncFileTypes;
	std::unordered_map<std::string, std::tuple<double,double>> implementToEffectiveLag;
	std::unordered_map<std::string, std::tuple<double,double>> effectiveToInactiveLag;
	std::unordered_map<std::string, std::tuple<double,double>> effectiveness;

	// control rules
	std::vector<controlRule> controlRules;
	//std::vector<flexibleControlRule> flexibleControlRules;
	std::vector<std::string> dcControlTypes;

	// reporting parameters
	std::tuple<double,double> indexReportLag;
	std::tuple<double,double> nonDCReportLag;
	std::tuple<double,double> dcReportLag;

	// Dangerous Contacts parameters
	bool dangerousContacts_on;
	std::map<DiseaseStatus, double> dcRiskScale; //Must be map and not unordered map to accomodate older compilers (older doesn't allow using enum as template arg for function std::hash which is used by unordered_map).
	double maxDCScale;

	// diagnostic parameters - maps keyed by diagnosticType name
	bool diagnostic_on;
	bool testSuspects_on;
	std::vector<std::string> diagnosticTypes;
	std::unordered_map<std::string, std::string> diagnosticScales;
	std::unordered_map<std::string, std::string> diagnosticConstraintFunctions;
	std::unordered_map<std::string, std::vector<double>> diagnosticConstraintFuncParams;
	std::unordered_map<std::string, std::vector<std::string>> diagnosticConstraintFuncFiles;
	std::unordered_map<std::string, std::vector<std::string>> diagnosticConstraintFuncFileTypes;
	std::unordered_map<std::string, std::tuple<double,double>> testStartToCompleteLag;
	std::unordered_map<std::string, std::tuple<double,double>> sensitivity;

	// diagnostic rules
	std::vector<diagnosticRule> diagnosticRules;
	std::vector<std::string> dcDiagnosticTypes;

	// suspecting parameters
	std::tuple<double,double> indexSuspectLag;
	std::tuple<double,double> nonDCSuspectLag;
	std::tuple<double,double> dcSuspectLag;

	//DIVA parameters
	//bool diva_on;
	//std::tuple<double,double> effVaxToDiagLag;

	//bTB slaughter surveillance
	std::string btb_slaughtershed_fname;
	std::string btb_slaughter_prop_fname;

	//BTB infection parameters [120 - 129]
	int btb_local_kernel_choice;
	int btb_wildlife_kernel_choice;
	std::vector<std::vector<double>> btb_local_kernel_parameters;
	std::vector<std::vector<double>> btb_wildlife_kernel_parameters;
	double btb_rel_wildlife_kernel_weight; //Must be between 0 and 1.
	double btb_mixed_kernel_scaling_factor = 1.0; //Scales up or down entire mixed kernel.
	std::string btb_wildlife_density_file;

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
		void readConfig(std::string&, bool);
};



inline const Parameters* File_manager::getParams()
{
	return &params;
}

#endif // File_manager_h
