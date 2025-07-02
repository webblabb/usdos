//© 2024 Colorado State University
/*
A county. Inherits from Region.
Call update_shipping_probabilities_USAMMv2(counties, k) to generate shipment probabilities,
where the argument counties is a vector containing pointers to all
counties (including self) and k is the shipment kernel object to use for
the probabilities. This function will populate an alias table within the
county, which is then used by the function get_shipment_destination_county() to
generate shipments. Each call to this function returns a destination county
for one shipment.

To create a complete county:
    (1) Construct with x,y and name.
    (2) Set the area of the county with set_area(double). This is important,
        otherwise the shipment probabilities will be wrong.
    (3) Call set_parent_state(State*) to set what state the county is in.
    (4) Add farms as pointers with add_premises or
        as a vector of pointers with set_farm.
    (5) When all counties have been created, call update_shipping_probabilities_USAMMv2 as above.
*/

#ifndef COUNTY_H
#define COUNTY_H

#include <string>
#include <unordered_map>
#include <map>
#include <set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Region.h"
#include "shared_functions.h"

extern int verboseLevel;

class Farm;
class Farm_type;
class Prem_class;
class Premsize_rep;
class State;
//class Shipment_kernel;
class USAMMv2_parameters;
class USAMMv3_parameters;

class County : public Region
{

public:
    County(std::string id, std::string kernel_str, int idx); //idx is the index of the county in the FIPS_vector.
    County(std::string id, double x, double y, std::string kernel_str);
    ~County();

    void update_shipping_probabilities_USAMMv2(std::vector<County*>& in_counties); //Incomplete
    void add_premises(Farm* in_farm, const std::vector<Farm_type*>& in_all_farm_types,
                      std::map<Farm_type*, std::set<Prem_class*>> in_all_prem_classes);
    void add_market(Farm_type*, double vol);
    void set_area(double in_area);
    void set_wildlife_density(double dens) { wildlife_density = dens; }
    void set_weights(std::vector<double> in_weights);
    void set_covariates_USAMMv2(std::map<Farm_type*, USAMMv2_parameters>& up_map);
    void set_covariates_USAMMv3(std::vector<USAMMv3_parameters>& up_vec);
    void update_covariate_weights_USAMMv2(std::map<Farm_type*, USAMMv2_parameters>& up_map,
                                          std::string time_period);
    void update_covariate_weights_USAMMv3(std::vector<USAMMv3_parameters>& up_vec,
                                          std::string time_period);

    void init_USAMMv3_shipment_vectors(std::vector<USAMMv3_parameters>& up_vec);
    void init_ft_pcl_vec(const std::map<Farm_type*, std::set<Prem_class*>>& ft_pcl_map);
    //Recalculates the premsize rep vector using the currently loaded USAMMv3 parameters.
    void update_cc_CVI_shipment_rate(std::vector<USAMMv3_parameters>& up_vec,
                                 Farm_type* ft, std::string time_period);

    void unset_shipping_probabilities(bool full_reset); //Sets the is_set_shipping flag to false, forcing a recalculation of shipping probabilities following an update of USAMM parameters.
    void normalize_shipping_weights(Farm_type* ft, double o_norm);
    void set_parent_state(State* target);
    void set_all_counties(std::vector<County*>* in_counties);
    void set_national_avg_farm_vol(Farm_type* ft, double vol);
    void set_national_avg_feedl_vol(Farm_type* ft, double vol);
    void set_national_avg_mkt_vol(Farm_type* ft, double vol);

    void set_outgoing_county_CVI_weight(size_t fty_idx, size_t o_pcl_idx, size_t d_pcl_idx, double v);

    double get_area(); //Inlined
    double get_wildlife_density() { return wildlife_density; }
    int get_fips_code() { return fips_code; }
    int get_idx() { return idx; }; //The index where this county is found in the FIPS_vector.
    size_t get_n_premises(); //Inlined
    size_t get_n_premises(Farm_type* ft); //Inlined
    int get_n_farms(Farm_type* ft);//Inlined, returns number of premises classed as ordinary farms of type ft - not feedlots, not markets. Just farms.
    int get_n_feedl(Farm_type* ft);//Inlined, returns number of premises classed as feedlots of type ft.
    int get_n_mkt();//Inlined, returns number of markets in county.
    double get_total_farm_vol(Farm_type* ft);
    double get_total_feedl_vol(Farm_type* ft);
    double get_total_mkt_vol();
    std::vector<Farm*>& get_premises(); //Inlined
    std::vector<Farm*>& get_premises(Farm_type* ft);
    const std::vector<Farm*>& get_premises_by_class(std::string tag) { return all_premises_by_class[tag]; };
    std::vector<Farm_type*> get_farm_types_present();
    bool has_premises_of_type_class(Farm_type* ft, Prem_class* pcl);
    std::vector<Prem_class*>& get_prem_classes_by_type_idx(int ft_idx) { return prem_classes_by_type[ft_idx]; }
    std::unordered_map<std::string, int> get_statuses(); //Inlined
    double get_county_ocov_weight(Farm_type* ft);
    double get_county_dcov_weight(Farm_type* ft);
    double get_o_market_weight(Farm_type* ft);
    double get_d_market_weight(Farm_type* ft);
    double get_o_unnormalized_prem_weight_sum(Farm_type* ft);
    double get_d_unnormalized_prem_weight_sum(Farm_type* ft);
    double get_outgoing_county_CVI_weight(size_t fty_idx, size_t o_pcl_idx, size_t d_pcl_idx);

    //Slaughter shipment destinations
    bool is_set_slaughter_probs() { return slaughter_probs_set; }
    void set_slaughter_probs(std::vector<double> probabilties, const std::vector<int>& facility_ids);
    int get_slaughter_destination();


    //USAMMv3 shipment generation
//    double get_d_CVI_shipment_rate_sum(Prem_class* d_pcl, std::string time_period,
//                                   USAMMv3_parameters& up);
//    double get_o_CVI_shipment_rate(Farm_type* ft, Prem_class* o_pcl,
//                               Prem_class* d_pcl, std::string time_period,
//                               std::vector<USAMMv3_parameters>& up_vec);
//    const std::vector<double>& get_d_CVI_shipment_rate_vec(Farm_type* ft, Prem_class* o_pcl,
//                                                       Prem_class* d_pcl, std::string time_period,
//                                                       std::vector<USAMMv3_parameters>& up_map);

    const std::vector<double>& get_eval_N_vec_origin_CVI(int fty_idx, int pcl_idx,
                                                         USAMMv3_parameters& up,
                                                         int time_period_idx);
    const std::vector<double>& get_eval_N_vec_dest_CVI(int fty_idx, int pcl_idx,
                                                       USAMMv3_parameters& up,
                                                       int time_period_idx);
    double get_within_county_adj_terms_CVI(int fty_idx, int pcl_idx,
                                           USAMMv3_parameters& up,
                                           int time_period_idx);

    const std::vector<double>& get_USAMMv3_origin_prem_weights_CVI(size_t fty_idx, Prem_class* o_pcl,
                                                               std::string time_period, USAMMv3_parameters& up);
    double get_USAMMv3_origin_prem_weight_sum_CVI(size_t fty_idx, Prem_class* o_pcl,
                                              std::string time_period, USAMMv3_parameters& up);
    const std::vector<double>& get_USAMMv3_dest_prem_weights_CVI(size_t fty_idx, Prem_class* d_pcl,
                                                             std::string time_period, USAMMv3_parameters& up);

    State* get_parent_state(); //Inlined
    County* get_shipment_destination_county(Farm_type* ft); //Generates a destination county for a shipment ORIGINATING from this county.
    Farm* get_shipment_destination_premises(Farm_type* ft); //Generates a destination premises for a shipment INBOUND to this county.

    double cov_weight_fun(std::vector<double> cov_values,
                          std::vector<double> cov_parameters);
    void calculate_centroid();
    void print_bools();
    bool is_initialized();

private:
    int verbose;

    gsl_rng* R;
    int fips_code = -1;
    std::string kernel_str;
    int idx = -1; //idx is the index of the county in the FIPS_vector.
    double area;
    double wildlife_density = 0.0;
    State* parent_state;
    std::vector<Farm*> member_premises;
    std::map<std::string, std::vector<Farm*>> all_premises_by_class;
    std::vector<Farm_type*> all_farm_types;
    size_t n_time_periods = 4;
    std::map<Farm_type*, std::vector<std::string>> county_ocov_names;
    std::map<Farm_type*, std::vector<double>> ocov_values;
    std::map<Farm_type*, std::vector<std::string>> county_dcov_names;
    std::map<Farm_type*, std::vector<double>> dcov_values;
    int tot_n_mkts;
    double tot_mkt_volume;
    std::map<Farm_type*, double> mkt_volumes;
    std::map<Farm_type*, double> county_ocov_weights;
    std::vector<double> county_dcov_weights;
    std::map<Farm_type*, double> o_unnormalized_prem_weight_sum;
    std::map<Farm_type*, double> d_unnormalized_prem_weight_sum;
    std::map<Farm_type*, double> weighted_avg_d_prem_weights;
    std::map<Farm_type*, std::vector<double>> receiver_weights; //Vectors of each farms' unnormalized weight, used to determine shipment receiver within county.
    std::map<Farm_type*, std::vector<Farm*>> receiver_weights_corr_farms; //Vectors of pointers to farms corresponding to above weights.
    std::map<Farm_type*, gsl_ran_discrete_t*> receiver_prob_distributions_by_type; //This is a gsl discrete probability distribution used to generate an outcome from all possible receivers.

    //The combined representation of all premises in the county when acting as receivers or senders.
    //v = the sum of all individual premises weights h^phi multiplied by state weight and county weight.
    //Used when generating shipments with USAMMv3 parameters only.
//    Vec_d_5d dest_premsize_rep_CVI_by_quarter; //[time_period_idx][Farm_type_idx][o_prem_class_idx][d_prem_clss_idx] = v
//    Vec_d_5d origin_premsize_rep_CVI_by_quarter; //[time_period_idx][Farm_type_idx][o_prem_class_idx][d_prem_clss_idx] = v

    //[Farm_type][sending_prem_class][rec_prem_class][receiving_county]
    //total rate for all pĶremises in this county of sending_prem_class sending to all
    //premises of rec_prem_class in rec_county. Includes everyting; inflow, outflow,
    //covariates, kernel, etc.
//    std::map<std::string, Vec_d_4d*> cc_CVI_shipment_rates_by_quarter;

    //[Farm_type][sending_prem_type][rec_prem_type] The rate of shipments originating
    //from this county to all other counties by prem_class combination of sending/rec premises.
    //The sum of the innermost vector in cc_shipment_rates. Used to see how many shipments
    //originate from this county every time-step in one go so that the entire county can
    //be skipped if there are no shipments. Updated together with cc_shipment_rates.
//    std::map<std::string, Vec_d_3d*> c_origin_CVI_shipment_rates_by_quarter;


    Vec_i_3d N_vec_CVI; //[fty][pcl][size_bin] Number of premises in each size bin and prem class
    Vec_d_4d weighted_o_N_vec_CVI; //[time_idx][fty][pcl][size_bin] N_vec_CVI multiplied by the weight (h^phi_O) for each corresponding bin size.
    Vec_d_4d weighted_d_N_vec_CVI;//[time_idx][fty][pcl][size_bin] N_vec_CVI multiplied by the weight (h^phi_D) for each corresponding bin size.
    Vec_d_3d within_county_adj_terms_CVI;//[time_idx][fty][pcl] The total weigt to be removed from the c_c sum of premsize rep when the county is both origin and destination to compensate for premises not sending to themselves.
    Vec_d_3d outgoing_county_CVI_weight; // [fty][opcl][dpcl] The sum of the weights to all other counties excluding this county's origin prem weighs. Used to calculate a single premises outgoing shipment rate (by multiplying this value with the premises individal size). Used for slaughter shipment rate in bTB slaughter surveillance.

    //The sum of the evaluated bin size destination weights multiplied by the corresponding
    //number of premises of that class in that bin size. [Farm_type idx][Prem_class idx] = v.
    Vec_d_2d d_CVI_lambda_sums;

    //The relative origin and destination weights among the premises within the county.
    //[Farm_type][pcl][prem_idx]
    Vec_d_3d internal_oprem_weights_CVI;
    Vec_d_3d internal_dprem_weights_CVI;
    Vec_d_2d internal_oprem_weight_sums_CVI;

    std::unordered_map<Farm_type*, std::vector<Farm*>> farms_by_type;
    std::unordered_map<Farm_type*, int> n_farms_by_type, n_feedl_by_type, n_mkt_by_type;
    std::vector<std::vector<size_t>> n_premises_of_type_class;
    std::vector<std::vector<Prem_class*>> prem_classes_by_type;
    std::unordered_map<Farm_type*, double> poisson_mean;
//    std::vector<Alias_table<County*>> shipping_probabilities; //By farm type index.
    std::vector<double*> shipping_probabilities;
    std::vector<size_t> n_outcomes;
    std::vector<County*>* shipping_outcomes;
    std::vector<County*>* all_counties;
    std::unordered_map<Farm_type*, double> national_avg_farm_vol, national_avg_feedl_vol, national_avg_mkt_vol;

    //Slaughter shipment generation
    bool slaughter_probs_set = false;
    gsl_ran_discrete_t* slaughter_facility_lookup_table = nullptr;
    std::vector<int> slaughter_facility_ids;

    bool county_initialized = false;
    bool is_set_area = false;
    bool is_set_state = false;
    bool is_set_shipment = false;

    void eval_N_vecs_for_fty(int fty_idx, int time_period_idx, USAMMv3_parameters& up,
                             Vec_d_2d& ft_weighted_o_N_vec, Vec_d_2d& ft_weighted_d_N_vec,
                             Vec_d_1d& ft_within_county_adj_terms);
    double pop_covariate(std::string cov_name, std::vector<std::string>& cov_names,
                         std::vector<double>& cov_values);

    virtual void set_initialized(bool& parameter);
    virtual void all_initialized();
};


inline double County::get_area()
{
    if(!is_set_area)
        not_initialized();

    return area;
}

inline size_t County::get_n_premises()
{
    return member_premises.size();
}

inline size_t County::get_n_premises(Farm_type* ft)
{
    return farms_by_type[ft].size();
}

inline std::vector<Farm*>& County::get_premises()
{
//    if(!county_initialized)
//        not_initialized();

    return member_premises;
}

inline State* County::get_parent_state()
{
//    if(!county_initialized)
//        not_initialized();

    return parent_state;
}

#endif // COUNTY_H
