#ifndef USAMMV3_PARAMETERS_H
#define USAMMV3_PARAMETERS_H

#include "Farm.h"

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class County;
class Farm_type;
class Parameters;

struct USAMMv3_parstruct
{
    public:
        USAMMv3_parstruct(Farm_type* fty, std::string period, int period_idx);
        ~USAMMv3_parstruct() {}

        Farm_type* fty;
        std::string period;
        int period_idx; //This is where this particular period is in the various USAMM_temporal vectors in the Parameters object (config file).
        std::map<std::pair<int, int>, double> c_map;
        std::map<std::pair<int, int>, double> kMu_map;
        std::map<std::pair<int, int>, double> kNu_map;
        std::unordered_map<int, double> phiO_map;
        std::unordered_map<int, double> phiD_map;
        std::unordered_map<std::string, double> dhalf_map;
        std::unordered_map<std::string, double> hratio_map;
        std::unordered_map<std::string, double> dhalfCOM_map;
        std::unordered_map<std::string, double> hratioCOM_map;
        std::unordered_map<std::string, double> outflow_map;
        std::unordered_map<std::string, double> inflow_map;
        std::unordered_map<std::string, double> outflowCOM_map;
        std::unordered_map<std::string, double> inflowCOM_map;
        std::unordered_map<std::string, double> ocov_map;
        std::unordered_map<std::string, double> dcov_map;
        std::unordered_map<std::string, double> a_map;
        std::unordered_map<std::string, double> b_map;
};

class USAMMv3_parameters
{
    typedef std::map<std::string, std::vector<double>> str_vec_map;
    typedef std::vector<double> Vec_d_1d;
    typedef std::vector<Vec_d_1d> Vec_d_2d;

    public:
        USAMMv3_parameters(const Parameters* parameters, Farm_type* fty,
                           std::set<Prem_class*> pclasses,
                           std::map<Prem_class*, double> avg_prem_sizes);
        ~USAMMv3_parameters();

        std::vector<Prem_class*> getPremClasses() { return prem_classes; }
        Farm_type* getFarmType() { return fty; }
        size_t get_n_time_periods() { return n_time_periods; }
        ///Returns the line of parameter values that was used in this USAMMv2_parameters object.
        std::string get_generation_string() { return generation_string; };

        ///Returns the names of the origin covariates that are loaded into this USAMMv2_parameters object.
        std::vector<std::string> get_county_ocov_names() { return ocov_par_names; }
        ///Returns the names of the destination covariates that are loaded into this USAMMv2_parameters object.
        std::vector<std::string> get_county_dcov_names() { return dcov_par_names; }
        ///Returns the COUNTY-specific origin covariates.
        std::vector<double> get_county_o_covs(County* c);
        ///Returns the COUNTY-specific destination covariates.
        std::vector<double> get_county_d_covs(County* c);

        double get_c(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_kMu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_kNu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_phi_O(Prem_class* pcl_O, std::string time_period);
        double get_phi_D(Prem_class* pcl_D, std::string time_period);
        std::pair<double, double> get_ab(std::string state, std::string time_period);
        double get_outflow(std::string state, std::string time_period);
        double get_inflow(std::string state, std::string time_period);
        double get_ocov(std::string ocov_name, std::string time_peroid);
        double get_dcov(std::string dcov_name, std::string time_peroid);
        double get_avg_prem_size(Prem_class* pcl) { return avg_prem_sizes.at(pcl); }
        std::vector<double> get_evaluated_origin_bin_weights(Prem_class* pcl, std::string period)
            { return origin_size_weights.at(period)[pcl->idx]; }
        std::vector<double> get_evaluated_destination_bin_weights(Prem_class* pcl, std::string period)
            { return destination_size_weights.at(period)[pcl->idx]; }
        std::map<Prem_class*, std::map<int, double>>& get_origin_binweight_lookup(std::string period)
            { return origin_size_weights_lookup.at(period); }
        std::map<Prem_class*, std::map<int, double>>& get_dest_binweight_lookup(std::string period)
            { return dest_size_weights_lookup.at(period); }

        static int farmSizeBinningFun(int p_original_size, Farm_type* ft, Prem_class* pcl);
        static std::set<int> getSizeBins(Farm_type* ft, Prem_class* pcl) { return USAMMv3_parameters::size_bins.at(ft).at(pcl); }
        static void clear_size_bins() { USAMMv3_parameters::size_bins.clear(); }

    private:
        const Parameters* config;
        Farm_type* fty;
        std::vector<Prem_class*> prem_classes;
        std::map<Prem_class*, double> avg_prem_sizes;
        size_t fty_idx;
        std::string fty_name;
        std::string post_fname;
        int n_metadata_lines;
        size_t last_sampled_line = 0;
        unsigned int n_lines_in_posterior_file = 0;
        std::map<std::string, int> pclass_tag_to_idx_map;
        std::unordered_map<std::string, USAMMv3_parstruct> parstruct_ptrs_by_period;
        std::vector<USAMMv3_parstruct*> partstr_ptrs; //Only used for access when debugging as Allinea-DDT won't let me look inside the above map for some reason.
        bool o_cov_loaded;
        bool d_cov_loaded;
        ///Stores the COUNTY-specific covariates. I.e number of something or slaughter volume etc.
        str_vec_map ocov_values;
        ///Stores the COUNTY-specific covariates. I.e number of something or slaughter volume etc.
        str_vec_map dcov_values;
        std::vector<std::string> ocov_par_names;
        std::vector<std::string> dcov_par_names;
        ///Stores the time periods found in the USAMMv3 posterior file (Q1, Q2, ... for instance).
        std::set<std::string> USAMM_time_periods;
        ///Stores the time periods found in the config file (Q1, Q2, ... for instance),
        ///these must match the ones the USAMMv3 posterior file.
        std::set<std::string> config_time_periods;
        size_t n_time_periods;
        ///Stores the line read from the posterior file.
        std::string generation_string;

        std::map<std::string, Vec_d_2d> origin_size_weights; //The size bins v evaluated as v^phi_o
        std::map<std::string, Vec_d_2d> destination_size_weights; //The size bins v evaluated as v^phi_d
        //A way to get the evaluated size of a premises of a particular prem_class and binned size.
        std::map<std::string, std::map<Prem_class*, std::map<int, double>>> origin_size_weights_lookup;
        std::map<std::string, std::map<Prem_class*, std::map<int, double>>> dest_size_weights_lookup;

        int sample_post(std::string fpath);
        std::pair<double, double> ab_fun(double dhalf, double hratio);
        ///Reads county-level covariates from a file. The covariate names
        ///in that file must match those in the parameter file.
        void initialize_covariates(std::string fname, str_vec_map& cov_map,
                                   std::vector<std::string>& header, bool& covariates_loaded);
        void evaluate_size_bins();
        size_t get_one_random_line(std::ifstream& f, std::string& res_line, bool header = true);
        bool temporal_name_exists(std::string temporal_name);

        static std::map<Farm_type*, std::map<Prem_class*, std::set<int>>> size_bins;
};

#endif // USAMMV3_PARAMETERS_H
