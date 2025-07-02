//© 2024 Colorado State University
#include "State.h"
#include "Farm.h"
#include "County.h"
#include "shared_functions.h"

#include <gsl/gsl_randist.h>

State::State(std::string id, int state_code, int usamm_version) :
    Region(id),
    state_code(state_code)
{
    type = "state";
    size_t seed = generate_distribution_seed();
    R = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(R, seed);
    this->setup_usamm_fun_pointers(usamm_version);
}

State::State(std::string id, double x, double y, int state_code, int usamm_version) :
    Region(id, x, y),
    state_code(state_code)
{
    type = "state";

    size_t seed = generate_distribution_seed();
    R = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(R, seed);
    this->setup_usamm_fun_pointers(usamm_version);
}

State::~State()
{
    gsl_rng_free(R);
}

void State::add_county(County* in_county)
{
    bool already_present = false;
    // can this be replaced with find? For readability - yes. For performance - doesn't matter.
    for(County* existing_county : member_counties)
    {
        if(existing_county == in_county)
        {
            already_present = true;
        }
    }

    if(!already_present)
    {
        member_counties.push_back(in_county);
        std::vector<Farm_type*> ftypes = in_county->get_farm_types_present();
        for(Farm_type* ft : ftypes)
        {
            farm_types_present.insert(ft);
        }
    }
}

void State::set_aCVI(double in_a, Farm_type* in_type)
{
    aCVI_map[in_type] = in_a;
}

void State::set_bCVI(double in_b, Farm_type* in_type)
{
    bCVI_map[in_type] = in_b;
}

void State::set_aCOM(double in_a, Farm_type* in_type)
{
    aCOM_map[in_type] = in_a;
}

void State::set_bCOM(double in_b, Farm_type* in_type)
{
    bCOM_map[in_type] = in_b;
}

void State::set_N(double in_N, Farm_type* in_type)
{
    N_map[in_type] = in_N;
}

void State::set_N_per_day(double in_daily_N, Farm_type* in_type)
{
    daily_N_map[in_type] = in_daily_N;
}

void State::set_N_rem(double in_rem_N, Farm_type* in_type)
{
    rem_N_map[in_type] = in_rem_N;
}

void State::set_null_lambda(double in_null_lambda, Farm_type* in_type)
{
    null_lambda_map[in_type] = in_null_lambda;
}

void State::set_s(double in_s, Farm_type* in_type)
{
    s_map[in_type] = in_s;
}

void State::set_std(double in_std, Farm_type* in_type)
{
    std_map[in_type] = in_std;
}

void State::set_kurt(double in_kurt, Farm_type* in_type)
{
    kurt_map[in_type] = in_kurt;
}

void State::set_inflowCVI(double in_inflow, Farm_type* in_type)
{
    size_t ft_idx = size_t(in_type->get_index());
    if(inflowCVI_vec.size() < ft_idx+1)
    {
        inflowCVI_vec.resize(ft_idx+1);
    }
    inflowCVI_vec[ft_idx] = in_inflow;
}

void State::set_outflowCVI(double in_outflow, Farm_type* in_type)
{
    outflowCVI_map[in_type] = in_outflow;
}

void State::set_inflowCOM(double in_inflow, Farm_type* in_type)
{
    size_t ft_idx = size_t(in_type->get_index());
    if(inflowCOM_vec.size() < ft_idx+1)
    {
        inflowCOM_vec.resize(ft_idx+1);
    }
    inflowCOM_vec[ft_idx] = in_inflow;
}

void State::set_outflowCOM(double in_outflow, Farm_type* in_type)
{
    outflowCOM_map[in_type] = in_outflow;
}

int State::get_code()
{
    return state_code;
}

std::string State::get_code_str()
{
    return std::to_string(this->get_code());
}

double State::get_aCVI(Farm_type* ft)
{
    return aCVI_map[ft];
}

double State::get_bCVI(Farm_type* ft)
{
    return bCVI_map[ft];
}

double State::get_aCOM(Farm_type* ft)
{
    return aCOM_map[ft];
}

double State::get_bCOM(Farm_type* ft)
{
    return bCOM_map[ft];
}

double State::get_N(Farm_type* ft)
{
    return N_map[ft];
}

double State::get_remaining_N(Farm_type* ft)
{
    return rem_N_map[ft];
}

std::unordered_map<Farm_type*, double> State::get_N_map()
{
    return N_map;
}

double State::get_null_lambda(Farm_type* ft)
{
    return null_lambda_map[ft];
}

double State::get_shipping_rate(Farm_type* ft)
{
    return shipping_rates.at(ft);
}

double State::get_s(Farm_type* ft)
{
    return s_map[ft];
}

double State::get_inflowCVI(Farm_type* ft)
{
    return inflowCVI_vec[ft->get_index()];
}

double State::get_outflowCVI(Farm_type* ft)
{
    return outflowCVI_map[ft];
}

double State::get_inflowCOM(Farm_type* ft)
{
    return inflowCOM_vec[ft->get_index()];
}

double State::get_outflowCOM(Farm_type* ft)
{
    return outflowCOM_map[ft];
}

std::unordered_map<Farm_type*, double> State::get_null_lambda_map()
{
    return null_lambda_map;
}

std::vector<County*> State::get_member_counties()
{
    return member_counties;
}

void State::set_shipping_parameters_USAMMv2(int USAMM_version, USAMMv2_parameters& usamm_par,
                                            Farm_type* ft, std::string time_period,
                                            size_t days_in_period, size_t days_rem)
{
    if(USAMM_version == 1)
    {
        this->set_shipping_parameters_v1(usamm_par, ft, time_period,
                                         days_in_period, days_rem);
    }
    else if(USAMM_version == 2)
    {
        this->set_shipping_parameters_v2(usamm_par, ft, time_period,
                                        days_rem);
    }
}

void State::set_shipping_parameters_USAMMv3(int USAMM_version, USAMMv3_parameters& usamm_par,
                                            Farm_type* ft, std::string time_period,
                                            size_t days_rem)
{
    this->set_shipping_parameters_v3(usamm_par, ft, time_period,
                                     days_rem);
}

void State::setup_usamm_fun_pointers(int usamm_version)
{
    if(usamm_version == 1)
    {
        n_ship_gen_fun = &State::generate_shipments_N;
    }
    else if(usamm_version == 2)
    {
        n_ship_gen_fun = &State::generate_shipments_rate;
    }
    else
    {
        n_ship_gen_fun = nullptr;
    }
}


void State::set_shipping_parameters_v1(USAMMv2_parameters& usamm_par,
                                       Farm_type* ft, std::string time_period,
                                       size_t days_in_period, size_t days_rem)
{
        double a = usamm_par.get_a(state_code, time_period);
        double b = usamm_par.get_b(state_code, time_period);
        double N = usamm_par.get_N(state_code, time_period);
        double daily_N = N / 365;
        double s = usamm_par.get_s(state_code, time_period);
        set_aCVI(a, ft);
        set_bCVI(b, ft);
        set_N(N, ft);
        set_N_per_day(daily_N, ft);
        set_s(s, ft);

        if(days_rem < days_in_period)
        {
            if(N_todo.find(ft) != N_todo.end() and
               N_todo[ft].find(time_period) != N_todo[ft].end())
            {
                //This period has already been visited and has something in the todo-map, use that as remaining shipments to be done.
                size_t rem_N = N_todo[ft][time_period];
                set_N_rem(double(rem_N), ft);
            }
            else
            {
                double daily_shipping_prob = double(days_rem) / double(days_in_period);
                unsigned int rem_N = gsl_ran_binomial(R, daily_shipping_prob, N);
                set_N_rem(double(rem_N), ft);
                N_todo[ft][time_period] = (unsigned int)(N+0.5) - rem_N; //This is what is left for the next time we enter this period.
            }
        }
        else
        {
            set_N_rem(N, ft);
        }
}

void State::set_shipping_parameters_v2(USAMMv2_parameters& usamm_par,
                                       Farm_type* ft, std::string time_period,
                                       size_t days_rem)
{
        double a = usamm_par.get_a(state_code, time_period);
        double b = usamm_par.get_b(state_code, time_period);
        double null_lambda = usamm_par.get_lambda(state_code, time_period);
        double s = usamm_par.get_s(state_code, time_period);
        double stdev = usamm_par.get_std(state_code, time_period);
        double kurt = usamm_par.get_kurt(state_code, time_period);
        set_aCVI(a, ft);
        set_bCVI(b, ft);
        set_null_lambda(null_lambda, ft);
        set_s(s, ft);
        set_std(stdev, ft);
        set_kurt(kurt, ft);
}

void State::set_shipping_parameters_v3(USAMMv3_parameters& usamm_par,
                                       Farm_type* ft, std::string time_period,
                                       size_t days_rem)
{
        std::pair<double, double> abCVI = usamm_par.get_abCVI(state_code, time_period);
        double inflowCVI = usamm_par.get_inflowCVI(state_code, time_period);
        double outlflowCVI = usamm_par.get_outflowCVI(state_code, time_period);
        set_aCVI(abCVI.first, ft);
        set_bCVI(abCVI.second, ft);
        set_inflowCVI(inflowCVI, ft);
        set_outflowCVI(outlflowCVI, ft);

        if(usamm_par.getUSAMMv3SpeciesCode() == "s") //Swine so add commuter parameters.
        {
            std::pair<double, double> abCOM = usamm_par.get_abCOM(state_code, time_period);
            double inflowCOM = usamm_par.get_inflowCOM(state_code, time_period);
            double outlflowCOM = usamm_par.get_outflowCOM(state_code, time_period);
            set_aCOM(abCOM.first, ft);
            set_bCOM(abCOM.second, ft);
            set_inflowCOM(inflowCOM, ft);
            set_outflowCOM(outlflowCOM, ft);
        }
}

int State::generate_shipments_N(Farm_type* ft, int days_rem)
{
    double todays_shipping_prob = 1.0 / double(days_rem);
    size_t n_shipments_rem = size_t(this->get_remaining_N(ft));
    unsigned int n_shipments_today = draw_binom(n_shipments_rem, todays_shipping_prob);
    this->set_N_rem(double(n_shipments_rem - n_shipments_today), ft);
    return n_shipments_today;
}

int State::generate_shipments_rate(Farm_type* ft, int days_rem)
{
    //Days rem is not used here, the shipping rate is already expressed per day.
    double daily_shipment_rate = this->get_shipping_rate(ft);
    return draw_poisson(daily_shipment_rate);
}


void State::reset_N_todo()
{
    N_todo.clear();
}

void State::normalize_shipping_weights(Farm_type* ft)
{
    double o_norm = 0.0;
    for(County* c : member_counties)
    {
        o_norm += c->get_o_unnormalized_prem_weight_sum(ft);
    }
    for(County* c : member_counties)
    {
        //Also sets each individual premises' (incl. markets) origin weight
        //to its unnormalized weight / the norm (sum across all premises and markets in state).
        c->normalize_shipping_weights(ft, o_norm);
    }
}

void State::update_shipping_rate(Farm_type* ft)
{
    double counties_oweight_sum = 0.0;
    for(County* c : member_counties)
    {
        //Removed * c->get_county_ocov_weight(ft) from this expr. since cov weight is included at
        //the farm level (see County::update_covariate_weights). I think this is correct.
        counties_oweight_sum += c->get_o_unnormalized_prem_weight_sum(ft);
    }
    shipping_rates[ft] = counties_oweight_sum * null_lambda_map.at(ft);
}

int State::generate_daily_shipments(Farm_type* ft, int days_rem)
{
    return (this->*n_ship_gen_fun)(ft, days_rem); //Call relevant function to generate daily shipments depending on usamm version.
}

void State::set_initialized(bool& parameter)
{
    parameter = true;
    all_initialized();
}

void State::all_initialized()
{
    if(is_set_id)
    {
        state_initialized = true;
    }
}

int State::get_n_premises() const
{
    int sum = 0;
    for (auto c : member_counties){
    	sum += c->get_n_premises();
    }
    return sum;
}

int State::get_n_premises(Farm_type* ft) const
{
    int sum = 0;
    for (auto c : member_counties){
    	sum += c->get_n_premises(ft);
    }
    return sum;
}

//int State::get_premises(std::vector<Farm*>& farm_v) const
//{
//    for(Farm_type* ft : farm_types_present)
//    {
//        get_premises(farm_v, ft);
//    }
//    return 0;
//}
//
//int State::get_premises(std::vector<Farm*>& farm_v, Farm_type* ft) const
//{
//    for(County* c : member_counties)
//    {
//        std::vector<Farm*> local_farms = c->get_premises(ft);
//        farm_v.insert(farm_v.end(), local_farms.begin(), local_farms.end());
//    }
//    return 0;
//}

double State::get_total_farm_weight(Farm_type* ft)
{
    double weight_sum = 0.0;
    for(County* c : member_counties)
    {
        weight_sum += c->get_o_unnormalized_prem_weight_sum(ft);
    }
    return weight_sum;
}

void State::print_bools()
{
	std::cout << "For " << id << std::endl; // Region member
	std::cout << "id = " << is_set_id << std::endl;
	std::cout << "position = " << is_set_position << std::endl;
	std::cout << "region = " << region_initialized << std::endl;
}
