#ifndef USAMM_PARAMETERS_H
#define USAMM_PARAMETERS_H

#include <vector>
#include <fstream>
#include <map>
#include <set>

struct Parameters;
class Farm_type;
class County;

//A map of maps for storing state.time_period.parameter_value.
typedef std::map<int, std::map<std::string, double>> int_par_map;
//Map of maps for storing time_period.covariate_name.parameter_value
typedef std::map<std::string, std::map<std::string, double>> str_par_map;
//Map of vectors for storing county.[covariate1, covariate2...]
typedef std::map<std::string, std::vector<double>> str_vec_map;

/*!
The USAMM_parameters class is used to read and store parameter output from
USAMM. The object is constructed with a pointer to a parameter struct where
information from the config file can be found. Also a farm type has to be
specified on construction. Provides functions to conveniently get parameter
values for different states and time periods.
This class also stores info about the county specific covariates (the actual
covariates, not only the covariate parameters) and provides access to those.
*/
class USAMM_parameters
{
    public:
        USAMM_parameters();
        USAMM_parameters(const Parameters* p, Farm_type* ft);
        ~USAMM_parameters();
        ///Reads parameter data from a USAMM .res file and stores them in a useful way.
        void initialize_parameters(std::string fname);
        ///Reads county-level covariates from a file. The covariate names
        ///in that file must match those in the parameter file.
        void initialize_covariates(std::string fname, str_vec_map& cov_map,
                                   std::vector<std::string>& header, bool& covariates_loaded);
        ///Reads a file containing info on what counties are considered super-nodes.
        void initialize_supernode_cov(std::string fname, str_vec_map& ocov_map,
                                      str_vec_map& dcov_map);
        ///Returns the farm type of this object.
        Farm_type* get_farm_type();
        double get_std(int state_code, std::string time_period);
        double get_kurt(int state_code, std::string time_period);
        double get_a(int state_code, std::string time_period);
        double get_b(int state_code, std::string time_period);
        double get_N(int state_code, std::string time_period);
        double get_lambda(int state_code, std::string time_period);
        double get_s(int state_code, std::string time_period);
        ///Returns the COUNTY-specific origin covariates.
        std::vector<double> get_county_o_covs(County* c);
        ///Returns the COUNTY-specific destination covariates.
        std::vector<double> get_county_d_covs(County* c);
        ///Returns the names of the origin covariates that are loaded into this USAMM_parameters object.
        std::vector<std::string> get_ocov_names();
        ///Returns the names of the destination covariates that are loaded into this USAMM_parameters object.
        std::vector<std::string> get_dcov_names();
        ///Returns the nation-wide origin covariate parameters.
        double get_ocov_parameter(std::string par_name, std::string time_period);
        ///Returns the nation-wide destination covariate parameters
        double get_dcov_parameter(std::string par_name, std::string time_period);
        ///Returns the line of parameter values that was used in this USAMM_parameters object.
        std::string get_generation_string();
    private:
        const Parameters* p;
        Farm_type* farm_type;
        size_t species_index;
        std::string species_name;
        bool supernodes_on;
        bool o_covariates_loaded;
        bool d_covariates_loaded;
        bool use_raw_parameters;
        ///Stores the COUNTY-specific covariates. I.e number of something or slaughter volume etc.
        str_vec_map county_ocov_map;
        ///Stores the COUNTY-specific covariates. I.e number of something or slaughter volume etc.
        str_vec_map county_dcov_map;
        std::vector<std::string> ocov_parameter_names;
        std::vector<std::string> dcov_parameter_names;
        int_par_map std_map;
        int_par_map kurt_map;
        int_par_map N_map;
        int_par_map lambda_map;
        int_par_map s_map;
        int_par_map a_map;
        int_par_map b_map;
        ///Stores the NATION-level PARAMETER values for different covariates, in other words the scaling factor for different covariates.
        str_par_map ocov_parameter_map;
        ///Stores the NATION-level PARAMETER values for different covariates, in other words the scaling factor for different covariates.
        str_par_map dcov_parameter_map;
        ///Stores the time periods found in the USAMM .res file (Q1, Q2, ... for instance).
        std::set<std::string> USAMM_time_periods;
        ///Stores the time periods found in the config file (Q1, Q2, ... for instance),
        ///these must match the ones the USAMM .res file.
        std::set<std::string> config_time_periods;
        ///Stores the line read from the .res file.
        std::string generation_string;

        ///Returns a and b (parameters used by the shipment distance kernel) from standard deviation
        ///and kurtosis as reported in the USAMM .res file.
        std::vector<double> calculate_a_and_b(double log_std, double log_kurt,
                                              double x1 = 0.5, double x2 = 0.05);
        ///Used for the calculation of a and b.
        bool get_parameters(double log_std, double log_kurt,double const & epsilon,
                            std::vector<double>& ab_vec);
        ///Used for the calculation of a and b.
        double log_int_hal(double n0,double n1, double epsilon, double log_kurtosis,
                           bool& status);
        ///Used for the calculation of a and b.
        double log_kurtosisf(double n, double log_kurtosis, bool& status);
        ///Used for the calculation of a and b.
        double return_log_a(double n, double sd);
        ///Used to select one random line (= set of parameters) from the USAMM output .res file.
        bool get_one_random_line(std::ifstream& f, std::string& res_line, bool header = true);
        bool temporal_name_exists(std::string temporal_name);
};

#endif // USAMM_PARAMETERS_H
