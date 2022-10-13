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

extern int verboseLevel;

class Farm;
class Farm_type;
class Prem_class;
class Premsize_rep;
class State;
class Shipment_kernel;
class USAMMv2_parameters;
class USAMMv3_parameters;

class County : public Region
{
    typedef std::vector<double> Vec_d_1d;
    typedef std::vector<Vec_d_1d> Vec_d_2d;
    typedef std::vector<Vec_d_2d> Vec_d_3d;
    typedef std::vector<Vec_d_3d> Vec_d_4d;

    typedef std::vector<int> Vec_i_1d;
    typedef std::vector<Vec_i_1d> Vec_i_2d;
    typedef std::vector<Vec_i_2d> Vec_i_3d;

public:
    County(std::string id, std::string kernel_str);
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
    void update_cc_shipment_rate(std::vector<USAMMv3_parameters>& up_vec,
                                 Farm_type* ft, std::string time_period);

    void unset_shipping_probabilities(bool full_reset); //Sets the is_set_shipping flag to false, forcing a recalculation of shipping probabilities following an update of USAMM parameters.
    void normalize_shipping_weights(Farm_type* ft, double o_norm);
    void set_parent_state(State* target);
    void set_all_counties(std::vector<County*>* in_counties);
    void set_national_avg_farm_vol(Farm_type* ft, double vol);
    void set_national_avg_feedl_vol(Farm_type* ft, double vol);
    void set_national_avg_mkt_vol(Farm_type* ft, double vol);

    double get_area(); //Inlined
    double get_wildlife_density() { return wildlife_density; }
    int get_fips_code() { return fips_code; }
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

    //Slaughter shipment destinations
    bool is_set_slaughter_probs() { return slaughter_probs_set; }
    void set_slaughter_probs(std::vector<double> probabilties, const std::vector<int>& facility_ids);
    int get_slaughter_destination();


    //USAMMv3 shipment generation
    double get_d_shipment_rate_sum(Prem_class* d_pcl, std::string time_period,
                                   USAMMv3_parameters& up);
    double get_o_shipment_rate(Farm_type* ft, Prem_class* o_pcl,
                               Prem_class* d_pcl, std::string time_period,
                               std::vector<USAMMv3_parameters>& up_vec);
    const std::vector<double>& get_d_shipment_rate_vec(Farm_type* ft, Prem_class* o_pcl,
                                                       Prem_class* d_pcl, std::string time_period,
                                                       std::vector<USAMMv3_parameters>& up_map);
    double get_rec_premsize_rep(Farm_type* ft, Prem_class* o_pcl, Prem_class* d_pcl,
                                std::vector<USAMMv3_parameters>& up_vec,
                                std::string time_period, int time_period_idx); //Returns the premsize representation of all premises of a certain Farm_type in the county. Used when calculating shipment rate from one arbitrary prem to all prems in this county. Includes prem size, phi and c.
    std::vector<double>& get_USAMMv3_origin_prem_weights(Farm_type* ft, Prem_class* o_pcl, std::string time_period,
                                                         std::vector<USAMMv3_parameters>& up_map);
    std::vector<double>& get_USAMMv3_dest_prem_weights(Farm_type* ft, Prem_class* d_pcl, std::string time_period,
                                                       std::vector<USAMMv3_parameters>& up_map);

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
    double area;
    double wildlife_density = 0.0;
    State* parent_state;
    std::vector<Farm*> member_premises;
    std::map<std::string, std::vector<Farm*>> all_premises_by_class;
    std::vector<Farm_type*> all_farm_types;
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

    //The combined representation of all premises in the county when acting as receivers
    //of a shipment from a premises of a specific prem class.
    //rec_premsize_rep[Farm_type_idx][sending_prem_class_idx] = v, where v = the sum
    // of all v^phiD * c * inflow * dcov_weights where c depends on the receiving premises class and the sending
    //premises class and phiD depends on the receiving premises class. Used when generating
    //shipments with USAMMv3 parameters only.
    Vec_d_4d rec_premsize_rep_by_quarter;

    //[Farm_type][sending_prem_class][rec_prem_class][receiving_county]
    //total rate for all premises in this county of sending_prem_class sending to all
    //premises of rec_prem_class in rec_county. Includes everyting; inflow, outflow,
    //covariates, kernel, etc.
    std::map<std::string, Vec_d_4d*> cc_shipment_rates_by_quarter;

    //[Farm_type][sending_prem_type][rec_prem_type] The rate of shipments originating
    //from this county to all other counties by prem_class combination of sending/rec premises.
    //The sum of the innermost vector in cc_shipment_rates. Used to see how many shipments
    //originate from this county every time-step in one go so that the entire county can
    //be skipped if there are no shipments. Updated together with cc_shipment_rates.
    std::map<std::string, Vec_d_3d*> c_origin_shipment_rates_by_quarter;

    //[Farm_type][prem_class][size_bin] Number of premises in each size bin and prem class
    Vec_i_3d N_vec;

    //The sum of the evaluated bin size destination weights multiplied by the corresponding
    //number of premises of that class in that bin size. [Farm_type idx][Prem_class idx] = v.
    Vec_d_2d d_lambda_sums;

    //The relative origin and destination weights among the premises within the county.
    //[Farm_type][pcl][prem_idx]
    std::map<Farm_type*, std::map<Prem_class*, std::vector<double>>> internal_oprem_weights;
    std::map<Farm_type*, std::map<Prem_class*, std::vector<double>>> internal_dprem_weights;

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
    if(!county_initialized)
        not_initialized();

    return member_premises;
}

inline State* County::get_parent_state()
{
    if(!county_initialized)
        not_initialized();

    return parent_state;
}

#endif // COUNTY_H
