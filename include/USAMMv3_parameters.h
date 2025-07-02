//© 2024 Colorado State University
#ifndef USAMMV3_PARAMETERS_H
#define USAMMV3_PARAMETERS_H

#include "Farm.h"
#include "shared_functions.h"

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

class County;
class Farm_type;
class Parameters;
class COM_prem_id;

struct USAMMv3_common_parstruct //Struct for parameters without temporal dimension.
{
    public:
        USAMMv3_common_parstruct() {};
        USAMMv3_common_parstruct(Farm_type* fty);
        ~USAMMv3_common_parstruct() {}

        Farm_type* fty;
        std::unordered_map<int, double> alpha_map; //Id and value. Ids are given by first line in id posterior file.
};

struct USAMMv3_period_parstruct //Struct for parameters with temporal dimension (i.e. quarters).
{
    public:
        USAMMv3_period_parstruct() {};
        USAMMv3_period_parstruct(Farm_type* fty, std::string period, int period_idx,
                                 USAMMv3_common_parstruct* common_pars);

        Farm_type* fty;
        std::string period;
        int period_idx; //This is where this particular period is in the various USAMM_temporal vectors in the Parameters object (config file).
        USAMMv3_common_parstruct* common_pars;
        std::vector<double> c_vec;
        std::vector<double> rho_vec;

        std::vector<double> kMu_vec;
        std::vector<double> kNu_vec;

        std::vector<double> shipSizeOScalingCVI_vec;
        std::vector<double> shipSizeDScalingCVI_vec;
        std::vector<double> shipSizeShapeCVI_vec;
        std::vector<double> shipSizeInterceptCVI_vec;

        std::vector<double> shipSizeOScalingCOM_vec;
        std::vector<double> shipSizeDScalingCOM_vec;
        std::vector<double> shipSizeShapeCOM_vec;
        std::vector<double> shipSizeInterceptCOM_vec;

        std::vector<double> phiO_vec;
        std::vector<double> phiD_vec;

        std::vector<double> dhalfCVI_vec;
        std::vector<double> hratioCVI_vec;
        std::vector<double> dhalfCOM_vec;
        std::vector<double> hratioCOM_vec;
        std::vector<double> outflowCVI_vec;
        std::vector<double> inflowCVI_vec;
        std::vector<double> outflowCOM_vec;
        std::vector<double> inflowCOM_vec;
        std::unordered_map<std::string, double> ocov_map;
        std::unordered_map<std::string, double> dcov_map;
        std::vector<double> aCVI_vec;
        std::vector<double> bCVI_vec;
        std::vector<double> aCOM_vec;
        std::vector<double> bCOM_vec;

        Vec_d_2d county_eval_kernel_CVI;
        Vec_d_2d county_eval_kernel_COM;
        Vec_d_4d eval_cc_rates_CVI; //Evaluated CVI shipment rates between all county pairs: [opcl][dpcl][ocounty][dcounty]
        Vec_d_4d eval_cc_rates_CVI_no_psize_rep; //Evaluated CVI shipment rates excluding the wight contributed by the premises themselves.
        Vec_d_3d eval_c_o_rates_CVI; //Evaluated totL CVI shipment rate out from a single county: [opcl][dpcl][ocounty]
        Vec_d_3d COM_random_effects; //Matrix of the log random effects for the rate between every pair of COM_ids. One for each farm type.
};

class USAMMv3_parameters
{
    public:
        USAMMv3_parameters(const Parameters* parameters, Farm_type* fty,
                           std::set<Prem_class*> pclasses,
                           std::map<Prem_class*, double> avg_prem_sizes,
                           std::vector<County*>& county_vec,
                           const Vec_d_2d& county_distance_matrix);
        ~USAMMv3_parameters();

        std::vector<Prem_class*> getPremClasses() { return prem_classes; }
        Farm_type* getFarmType() { return fty; }
        std::string getUSAMMv3SpeciesCode() { return USAMMv3_species_code; }
        std::map<int, int> get_unobserved_premises() {return unobserved_prem_county_association; }
        size_t get_n_time_periods() { return n_time_periods; }
        std::string get_generation_string() { return generation_string; };
        USAMMv3_period_parstruct* get_parstruct_by_period(std::string time_period) { return parstructs_by_period.at(time_period); }
        USAMMv3_period_parstruct& get_parstruct_by_period_idx(int time_period_idx) { return parstruct_vec.at(time_period_idx); }

        ///Returns the names of the origin covariates that are loaded into this USAMMv2_parameters object.
        std::vector<std::string> get_county_ocov_names() { return ocov_par_names; }
        ///Returns the names of the destination covariates that are loaded into this USAMMv2_parameters object.
        std::vector<std::string> get_county_dcov_names() { return dcov_par_names; }
        ///Returns the COUNTY-specific origin covariates.
        std::vector<double> get_county_o_covs(County* c);
        ///Returns the COUNTY-specific destination covariates.
        std::vector<double> get_county_d_covs(County* c);

        double get_c(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_rho(std::string time_period);

        double get_kMu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_kNu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);

        double get_oShipSizeScalingCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_dShipSizeScalingCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_shipSizeShapeCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_shipSizeInterceptCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_oShipSizeScalingCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_dShipSizeScalingCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_shipSizeShapeCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);
        double get_shipSizeInterceptCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period);

        double get_phi_O(Prem_class* pcl_O, std::string time_period);
        double get_phi_D(Prem_class* pcl_D, std::string time_period);
        std::pair<double, double> get_abCVI(int, std::string time_period);
        std::pair<double, double> get_abCOM(int, std::string time_period);
        double get_outflowCVI(int, std::string time_period);
        double get_outflowCOM(int, std::string time_period);
        double get_inflowCVI(int, std::string time_period);
        double get_inflowCOM(int, std::string time_period);
        double get_ocov(std::string ocov_name, std::string time_peroid);
        double get_dcov(std::string dcov_name, std::string time_peroid);
        double get_producer_alpha(int producer_id);
        std::vector<COM_prem_id>& get_COM_ids() { return COM_ids; }
        std::vector<COM_prem_id*>& get_COM_id_ptrs_by_producer(int producer_idx) { return COM_ids_by_producer_idx.at(producer_idx); }
        double get_avg_prem_size(Prem_class* pcl) { return avg_prem_sizes.at(pcl); }
        const Vec_d_2d& get_county_distance_matrix() { return county_distance_matrix; }
        double kernel_fun(double d, double a, double b);
        std::vector<double> get_evaluated_origin_bin_weights(Prem_class* pcl, int period_idx)
            { return origin_size_weights.at(period_idx)[pcl->idx]; }
        std::vector<double> get_evaluated_destination_bin_weights(Prem_class* pcl, int period_idx)
            { return destination_size_weights.at(period_idx)[pcl->idx]; }
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
        std::string USAMMv3_species_code;
        std::map<Prem_class*, double> avg_prem_sizes;
        std::vector<County*>& counties;
        const Vec_d_2d& county_distance_matrix;
        size_t fty_idx;
        std::string fty_name;
        std::string post_fname;
        std::string COM_id_fpath = "";
        int n_metadata_lines;
        size_t last_sampled_line = 0;
        unsigned int n_lines_in_posterior_file = 0;
        std::map<std::string, int> pclass_tag_to_idx_map;
        USAMMv3_common_parstruct common_parstruct;
        std::unordered_map<std::string, USAMMv3_period_parstruct*> parstructs_by_period;
        std::vector<USAMMv3_period_parstruct> parstruct_vec;
        std::vector<USAMMv3_period_parstruct*> partstr_ptrs; //Only used for access when debugging as Allinea-DDT won't let me look inside the above map for some reason.
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
        std::vector<std::string> USAMM_time_period_vec;
        ///Stores the time periods found in the config file (Q1, Q2, ... for instance),
        ///these must match the ones the USAMMv3 posterior file.
        std::set<std::string> config_time_periods;
        size_t n_time_periods;
        ///Stores the line read from the posterior file.
        std::string generation_string;

        std::map<std::string, int> COM_producers_name_to_id;
        std::vector<COM_prem_id> COM_ids;
        std::vector<std::vector<COM_prem_id*>> COM_ids_by_producer_idx;
        std::map<int, COM_prem_id*> COM_ids_by_id;
        std::map<int, int> unobserved_prem_county_association;
        std::vector<int> unobserved_prem_ids;
        std::map<int, int> unobs_prem_sizes;

        std::map<int, Vec_d_2d> origin_size_weights; //The size bins v evaluated as v^phi_o
        std::map<int, Vec_d_2d> destination_size_weights; //The size bins v evaluated as v^phi_d
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
        void split_pclasses(const std::string& parname, int& pcl_idx_1, int& pcl_idx_2);
        int readSwineExclusiveData(std::string fpath, size_t posterior_line);

        static std::map<Farm_type*, std::map<Prem_class*, std::set<int>>> size_bins;
        std::unordered_map<std::string, int> abbr_to_fips_state_map =
            { {"AL", 1}, {"AZ", 4}, {"AR", 5}, {"CA", 6}, {"CO", 8}, {"CT", 9}, {"DE", 10}, {"FL", 12},
              {"GA", 13}, {"ID", 16}, {"IL", 17}, {"IN", 18}, {"IA", 19}, {"KS", 20}, {"KY", 21},
              {"LA", 22}, {"ME", 23}, {"MD", 24}, {"MA", 25}, {"MI", 26}, {"MN", 27}, {"MS", 28},
              {"MO", 29}, {"MT", 30}, {"NE", 31}, {"NV", 32}, {"NH", 33}, {"NJ", 34}, {"NM", 35},
              {"NY", 36}, {"NC", 37}, {"ND", 38}, {"OH", 39}, {"OK", 40}, {"OR", 41}, {"PA", 42},
              {"RI", 44}, {"SC", 45}, {"SD", 46}, {"TN", 47}, {"TX", 48}, {"UT", 49}, {"VT", 50},
              {"VA", 51}, {"WA", 53}, {"WV", 54}, {"WI", 55}, {"WY", 56} };
};

class COM_prem_id
{
public:
    COM_prem_id(int id);
    ~COM_prem_id();
//    void initialize_random_eff_vectors(size_t n_com_ids, size_t n_periods);
    int get_id() { return id; }
    std::string get_id_str() { return id_str; }
    size_t get_idx() { return idx; }
    int get_sampled_prem_id() { return sampled_prem_id; }
    Farm* get_farm() { return farm; }
    County* get_county() { return county; }
    int get_n_shipments_to(COM_prem_id* other, size_t period_idx);
    double get_random_effect(size_t other_COM_id_idx, size_t period_idx);
    void update_farm(Farm* in_farm);
    void set_sampled_prem_id(int id) { sampled_prem_id = id; } //Sets the id-number for the prem that his COM_id is to be associated to.
    void set_idx(size_t p_idx) { idx = p_idx; }
    void set_n_shipments_to(COM_prem_id* other, int period_idx, int n) { n_shipments_to_others_by_period.at(period_idx)[other] = n; }
//    void set_random_effect(size_t other_COM_id_idx, size_t period_idx, double val);
    void set_county(County* in_county) { county = in_county; }
    void set_producer(int p);
    int get_producer() { return com_producer; }
    void disconnect_farm();


private:
    int id;
    size_t idx = 999999999;
    std::string id_str;
    int sampled_prem_id = 999999999;
    Farm* farm = nullptr;
    County* county = nullptr;
    std::vector<std::map<COM_prem_id*, int>> n_shipments_to_others_by_period; //Number of shipments to other COM_prem_ids (only the ones with > 0 are in the map).

    int com_producer = -1;
};

#endif // USAMMV3_PARAMETERS_H
