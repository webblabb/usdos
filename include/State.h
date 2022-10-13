/* A state, inherits from region */

#ifndef STATE_H
#define STATE_H

#include <string>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include "Region.h"
#include "USAMMv2_parameters.h"
#include "USAMMv3_parameters.h"
#include <gsl/gsl_rng.h>

class Farm;
class County;
class Farm_type;

class State : public Region
{
    typedef int (State::*n_ship_gen_fun_ptr)(Farm_type*, int);

public:
    State(std::string id, int state_code, int usamm_version = 0);
    State(std::string name, double x, double y, int state_code, int usamm_version = 0);
    ~State();

    void add_county(County* in_county);
    void set_a(double in_a, Farm_type* in_type);
    void set_b(double in_b, Farm_type* in_type);
    void set_N(double in_N, Farm_type* in_type);
    void set_N_per_day(double in_daily_N, Farm_type* in_type);
    void set_N_rem(double in_rem_N, Farm_type* in_type);
    void set_null_lambda(double in_null_lambda, Farm_type* in_type);
    void set_s(double in_s, Farm_type* in_type);
    void set_std(double in_std, Farm_type* in_type);
    void set_kurt(double in_kurt, Farm_type* in_type);
    void set_inflow(double in_inflow, Farm_type* in_type);
    void set_outflow(double in_outflow, Farm_type* in_type);
    void set_shipping_parameters_USAMMv2(int USAMM_version, USAMMv2_parameters& usamm_par,
                                         Farm_type* ft, std::string time_period,
                                         size_t days_in_period, size_t days_rem);
    void set_shipping_parameters_USAMMv3(int USAMM_version, USAMMv3_parameters& usamm_par,
                                         Farm_type* ft, std::string time_period,
                                         size_t days_rem);
    void reset_N_todo();
    void normalize_shipping_weights(Farm_type* ft);
    void update_shipping_rate(Farm_type* ft);
    int generate_daily_shipments(Farm_type* ft, int days_rem);

    int get_code();
    std::string get_code_str();
    double get_a(Farm_type* ft);
    double get_b(Farm_type* ft);
    double get_N(Farm_type* ft);
    double get_remaining_N(Farm_type* ft);
    std::unordered_map<Farm_type*, double> get_N_map();
    double get_null_lambda(Farm_type* ft);
    double get_shipping_rate(Farm_type* ft); //This returns the sum of all member farms' weights multiplied by null_lambda.
    double get_s(Farm_type* ft);
    double get_inflow(Farm_type* ft);
    double get_outflow(Farm_type* ft);
    double get_total_farm_oweight_sum(); //Needed when figuring out probability that shipments originate from unaffected farms ( P(unaffected) = P(Total weight) - P(affected) ).
    std::unordered_map<Farm_type*, double> get_null_lambda_map();
    std::vector<County*> get_member_counties();
    int get_n_counties(); //inlined
    int get_n_premises() const;
    int get_n_premises(Farm_type* ft) const;
//    int get_premises(std::vector<Farm*>& farm_v) const;
//    int get_premises(std::vector<Farm*>& farm_v, Farm_type* ft) const;
    double get_total_farm_weight(Farm_type* ft);
    void print_bools();

private:
    int state_code;
    gsl_rng* R;
    std::vector<County*> member_counties;
    std::unordered_set<Farm_type*> farm_types_present;
    bool state_initialized = false;
    std::unordered_map<Farm_type*, double> a_map;
    std::unordered_map<Farm_type*, double> b_map;
    std::unordered_map<Farm_type*, double> N_map;
    std::unordered_map<Farm_type*, double> daily_N_map;
    //rem_N_map stores what is left between the current timestep and the end of the current time period (e.g. quarter).
    std::unordered_map<Farm_type*, double> rem_N_map;
    //N_todo stores what remaing to be done for the time period. If the simulation/generation starts in the middle
    //of a period, there will remain some amount of shipments to be done when that period is entered the next time.
    std::unordered_map<Farm_type*, std::unordered_map<std::string, size_t>> N_todo;
    std::unordered_map<Farm_type*, double> null_lambda_map;
    std::unordered_map<Farm_type*, double> shipping_rates;
    std::unordered_map<Farm_type*, double> s_map;
    std::unordered_map<Farm_type*, double> std_map;
    std::unordered_map<Farm_type*, double> kurt_map;
//    std::unordered_map<Farm_type*, double> inflow_map;
    std::vector<double> inflow_vec;
    std::unordered_map<Farm_type*, double> outflow_map;

    virtual void set_initialized(bool& parameter);
    virtual void all_initialized();

    void setup_usamm_fun_pointers(int usamm_version);
    void set_shipping_parameters_v1(USAMMv2_parameters& usamm_par, Farm_type* ft, //For use with usamm v1 parameters
                                    std::string time_period, size_t days_in_period, size_t days_rem);
    void set_shipping_parameters_v2(USAMMv2_parameters& usamm_par, Farm_type* ft, //For use with usamm v2 parameters
                                    std::string time_period, size_t days_rem);
    void set_shipping_parameters_v3(USAMMv3_parameters& usamm_par, Farm_type* ft, //For use with usamm v2 parameters
                                    std::string time_period, size_t days_rem);

    int generate_shipments_N(Farm_type* ft, int days_rem);
    int generate_shipments_rate(Farm_type* ft, int days_rem);

    n_ship_gen_fun_ptr n_ship_gen_fun; //Pointer to function that generates number of shipments for given day based on USAMM version (using N or lambda).
};

inline int State::get_n_counties()
{
    return int(member_counties.size());
}

#endif // STATE_H
