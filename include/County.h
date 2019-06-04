/*
A county. Inherits from Region.
Call update_shipping_probabilities(counties, k) to generate shipment probabilities,
where the argument counties is a vector containing pointers to all
counties (including self) and k is the shipment kernel object to use for
the probabilities. This function will populate an alias table within the
county, which is then used by the function get_shipment_destination() to
generate shipments. Each call to this function returns a destination county
for one shipment.

To create a complete county:
    (1) Construct with x,y and name.
    (2) Set the area of the county with set_area(double). This is important,
        otherwise the shipment probabilities will be wrong.
    (3) Call set_parent_state(State*) to set what state the county is in.
    (4) Add farms as pointers with add_farm or
        as a vector of pointers with set_farm.
    (5) When all counties have been created, call update_shipping_probabilities as above.
*/

#ifndef COUNTY_H
#define COUNTY_H

#include <string>
#include <unordered_map>
#include <map>
#include <gsl_rng.h>
#include "Region.h"
#include "Alias_table.h"

extern int verboseLevel;

class Farm;
class Farm_type;
class State;
class Shipment_kernel;
class USAMM_parameters;

class County : public Region
{
public:
    County(std::string id, std::string kernel_str);
    County(std::string id, double x, double y, std::string kernel_str);
    ~County();

    void set_farms(const std::vector<Farm*>& in_farms);
    void add_farm(Farm* in_farm);
    void update_shipping_probabilities(std::vector<County*>& in_counties); //Incomplete
    void set_area(double in_area);
    void set_weights(std::vector<double> in_weights);
    void set_covariates(std::map<Farm_type*, USAMM_parameters>& up_map);
    void update_covariate_weights(std::map<Farm_type*, USAMM_parameters>& up_map,
                                  std::string time_period);
    void unset_shipping_probabilities(); //Sets the is_set_shipping flag to false, forcing a recalculation of shipping probabilities follwing an update of USAMM parameters.
    void normalize_shipping_weight(Farm_type* ft, double norm);
    void set_parent_state(State* target);
    void set_all_counties(std::vector<County*>* in_counties);

    double get_area(); //Inlined
    size_t get_n_farms(); //Inlined
    size_t get_n_farms(Farm_type* ft); //Inlined
    std::vector<Farm*> get_farms(); //Inlined
    std::vector<Farm*> get_farms(Farm_type* ft);
    std::vector<Farm_type*> get_farm_types_present();
    std::unordered_map<std::string, int> get_statuses(); //Inlined
    double get_ocov_weight(Farm_type* ft);
    double get_dcov_weight(Farm_type* ft);
    double get_unnormed_o_weight(Farm_type* ft);
    double get_normed_o_weight(Farm_type* ft);
    double get_o_farm_weight_sum(Farm_type* ft);
    double get_d_farm_weight_sum(Farm_type* ft);

    State* get_parent_state(); //Inlined
    County* get_shipment_destination(Farm_type* ft);

    double cov_weight_fun(std::vector<double> cov_values,
                          std::vector<double> cov_parameters);
    void calculate_centroid();
    void print_bools();
    bool is_initialized();

private:
    int verbose;

    gsl_rng* R;
    std::string kernel_str;
    double area;
    State* parent_state;
    std::vector<Farm*> member_farms;
    std::mt19937 mt19937_generator;
    std::map<Farm_type*, std::vector<std::string>> ocov_names;
    std::map<Farm_type*, std::vector<double>> ocov_values;
    std::map<Farm_type*, std::vector<std::string>> dcov_names;
    std::map<Farm_type*, std::vector<double>> dcov_values;
    std::map<Farm_type*, double> ocov_weights;
    std::map<Farm_type*, double> dcov_weights;
    std::map<Farm_type*, double> normed_o_weight;
    std::map<Farm_type*, double> o_farm_weight_norms;
    std::map<Farm_type*, double> d_farm_weight_norms;
    std::map<Farm_type*, double> weighted_avg_d_farm_weights;
    std::unordered_map<Farm_type*, std::vector<Farm*>> farms_by_type;
    std::unordered_map<Farm_type*, double> poisson_mean;
    std::unordered_map<Farm_type*, std::poisson_distribution<int>*> poisson_map;
//    std::vector<Alias_table<County*>> shipping_probabilities; //By farm type index.
    std::vector<double*> shipping_probabilities;
    std::vector<size_t> n_outcomes;
    std::vector<County*>* shipping_outcomes;
    std::vector<County*>* all_counties;

    bool county_initialized = false;
    bool is_set_area = false;
    bool is_set_state = false;
    bool is_set_shipment = false;

    virtual void set_initialized(bool& parameter);
    virtual void all_initialized();
};

inline double County::get_area()
{
    if(!is_set_area)
        not_initialized();

    return area;
}

inline size_t County::get_n_farms()
{
    return member_farms.size();
}

inline size_t County::get_n_farms(Farm_type* ft)
{
    return farms_by_type[ft].size();
}

inline std::vector<Farm*> County::get_farms()
{
    if(!county_initialized)
        not_initialized();

    return member_farms;
}

inline State* County::get_parent_state()
{
    if(!county_initialized)
        not_initialized();

    return parent_state;
}

#endif // COUNTY_H
