#include "County.h"
#include "Farm.h"
#include "State.h"
#include "Shipment_kernel.h"
#include "USAMM_parameters.h"
#include "shared_functions.h"

#include <limits>
#include <algorithm>
#include <gsl_randist.h>

County::County(std::string id, std::string kernel_str) :
    Region(id), kernel_str(kernel_str), area(0.0)
//    mt19937_generator(generate_distribution_seed())
{
    verbose = verboseLevel;
    type = "county";
    size_t seed = generate_distribution_seed();
    R = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(R, seed);
}

County::County(std::string id, double x, double y, std::string kernel_str) :
    Region(id, x, y), kernel_str(kernel_str), area(0.0)
//    mt19937_generator(generate_distribution_seed())
{
    type = "county";
}

County::~County()
{
    for(size_t i = 0; i < shipping_probabilities.size(); i++)
    {
        delete shipping_probabilities.at(i);
    }
    gsl_rng_free(R);
}

//Measures the distance to all counties and calculates probabilities to send to them.
//Arguments: a vector of all th counties, a pointer to a function describing the
//kernel and a function that calculates the distance between two point objects (pointers).
void County::update_shipping_probabilities(std::vector<County*>& in_counties)
{
    for(size_t i = 0; i < shipping_probabilities.size(); i++)
    {
        delete shipping_probabilities.at(i);
    }

    if(!is_initialized())
    {
        std::cout << "Make sure all of the following are set before attempting to "
                  << "initialize shipment probabilities." << std::endl;
        print_bools();
        not_initialized();
    }

    shipping_probabilities.resize(farms_by_type.size());
    shipping_outcomes = &in_counties;
    n_outcomes.resize(farms_by_type.size());
    //Calculate shipping probabilities for all farm-types (species).
    for(auto ft_vec_pair : farms_by_type)
    {
        Farm_type* current_ft = ft_vec_pair.first;
        size_t ft_i = current_ft->Farm_type::get_index();
        size_t n = in_counties.size();
        n_outcomes.at(ft_i) = n;
        shipping_probabilities.at(ft_i) = new double[n];

        std::vector<double> probabilities;
        probabilities.reserve(in_counties.size());
        //The kernel is based on the states current a and b and gives the shape of the distance dependence.
        Shipment_kernel k(this->get_parent_state()->get_a(current_ft),
                          this->get_parent_state()->get_b(current_ft),
                          kernel_str, true);
        double normalization_sum = 0.0;

        //Get all kernel values and keep track of the total for use when normalizing.
        for(auto c : in_counties)
        {
            //Destination inflow weight based on number of possible farms that can receive in the receiving county,
            //therefore n-1 if sending within the county itself.
            double d_farm_weight = c->get_d_farm_weight_sum(current_ft);
            if(c == this)
            {
                d_farm_weight -= this->weighted_avg_d_farm_weights.at(current_ft);
            }

            double inflow_weight = c->get_parent_state()->get_s(current_ft) * d_farm_weight;
            //Destination covariate weight.
            double dcov_weight = c->get_dcov_weight(current_ft);
            //Kernel value * Flow of state of 'origin' county * number of farms in 'origin' county.
            double kernel_value = k.kernel(this, c);
            double unnormalized_probability = kernel_value * inflow_weight * dcov_weight;
            if(unnormalized_probability == 0.0 and
               inflow_weight != 0.0 and
               dcov_weight != 0.0)
            {
                unnormalized_probability = std::numeric_limits<double>::min();
            }

            probabilities.push_back(unnormalized_probability);
            normalization_sum += unnormalized_probability;
        }
        if(normalization_sum == 0.0)
        {
            std::cout << "The shipping probabilities of county " << this->get_id() <<
                         " are all zero. Exiting." << std::endl;
            exit(1);
        }

        for(size_t i = 0; i < probabilities.size(); i++)
        {
            shipping_probabilities.at(ft_i)[i] = probabilities[i] / normalization_sum;
        }

    }
    set_initialized(is_set_shipment);
}

//Sets the farms that belong to this county by passing a
//vector of pointers to them
void County::set_farms(const std::vector<Farm*>& in_farms)
{
    for(Farm* in_farm : in_farms)
    {
        this->add_farm(in_farm);
    }
}

//Adds one single farm that belongs to this county by passing a pointer to it.
//Also called internally by set_farms
void County::add_farm(Farm* in_farm)
{
    member_farms.push_back(in_farm);
    farms_by_type[in_farm->get_farm_type()].push_back(in_farm);
    in_farm->set_parent_county(this);
}

void County::set_area(double in_area)
{
    area = in_area;
    set_initialized(is_set_area);
}

void County::set_covariates(std::map<Farm_type*, USAMM_parameters>& up_map)
{
    for(auto& name_up_pair : up_map)
    {
        Farm_type* up_ft = name_up_pair.second.get_farm_type();
        ocov_names[up_ft] = name_up_pair.second.get_ocov_names();
        ocov_values[up_ft] = name_up_pair.second.get_county_o_covs(this);
        dcov_names[up_ft] = name_up_pair.second.get_dcov_names();
        dcov_values[up_ft] = name_up_pair.second.get_county_d_covs(this);
    }
}

void County::update_covariate_weights(std::map<Farm_type*, USAMM_parameters>& up_map,
                                      std::string time_period)
{
    for(auto& name_up_pair : up_map)
    {
        USAMM_parameters& up = name_up_pair.second;
        Farm_type* up_ft = up.get_farm_type();

        //Only calculate weights if there are any farms of the current farm type.
        if(this->get_n_farms(up_ft) > 0)
        {
            //Origin county-level covariates.
            std::vector<double> ocov_parameters(ocov_names.at(up_ft).size(), 0.0);
            for(size_t i = 0; i< ocov_parameters.size(); i++)
            {
                ocov_parameters[i] = up.get_ocov_parameter(ocov_names.at(up_ft).at(i), time_period);
            }
            ocov_weights[up_ft] = cov_weight_fun(ocov_values.at(up_ft), ocov_parameters);

            //Destination county-level covariates.
            std::vector<double> dcov_parameters(dcov_names.at(up_ft).size(), 0.0);
            for(size_t i = 0; i< dcov_parameters.size(); i++)
            {
                dcov_parameters[i] = up.get_dcov_parameter(dcov_names.at(up_ft).at(i), time_period);
            }
            dcov_weights[up_ft] = cov_weight_fun(dcov_values.at(up_ft), dcov_parameters);

            //Origin and destination farm level weight together.
            double temp_oweight_sum = 0.0;
            double temp_dweight_sum = 0.0;
            double o_farm_size_covariate = 0.0;//These are temporary placeholders, to be exchanged for
            double d_farm_size_covariate = 0.0;//real farm-level covariates at a later point.
            //Calculate the farm weight norms (sums) and set farm's weights.
            double temp_avg_farm_weight = 0.0;
            for(Farm* f : this->farms_by_type.at(up_ft))
            {
                double n_animals = double(f->get_size(up_ft->get_species()));
                double ow = 0.0;
                double dw = 0.0;
                if(n_animals > 0.0)
                {
                    ow = std::pow(n_animals, o_farm_size_covariate); //This is the farm-level covariate function. Maybe break it out as its own?
                    ow = ow * ocov_weights[up_ft]; //Weigh individual farm weight by county-level origin covariate weight.
                    dw = std::pow(n_animals, d_farm_size_covariate);
                }
                f->set_unnormalized_oweight(ow);
                f->set_unnormalized_dweight(dw);
                temp_oweight_sum += ow;
                temp_dweight_sum += dw;
                temp_avg_farm_weight += ow*dw;
            }
            o_farm_weight_norms[up_ft] = temp_oweight_sum;
            d_farm_weight_norms[up_ft] = temp_dweight_sum;
            weighted_avg_d_farm_weights[up_ft] = temp_avg_farm_weight / temp_oweight_sum;
        }
        //If not, set weights to 0.
        else
        {
            ocov_weights[up_ft] = 0.0;
            dcov_weights[up_ft] = 0.0;
            o_farm_weight_norms[up_ft] = 0.0;
            d_farm_weight_norms[up_ft] = 0.0;
            weighted_avg_d_farm_weights[up_ft] = 0.0;
        }
    }
}

void County::unset_shipping_probabilities()
{
    is_set_shipment = false;
}

void County::normalize_shipping_weight(Farm_type* ft, double norm)
{
    normed_o_weight[ft] = this->get_unnormed_o_weight(ft) / norm;
    for(Farm* f : farms_by_type.at(ft))
    {
        f->set_normalized_oweight(f->get_unnormalized_oweight() / norm);
    }
}

void County::set_parent_state(State* target)
{
    parent_state = target;
    target->add_county(this);
    set_initialized(is_set_state);
}

void County::set_all_counties(std::vector<County*>* in_counties)
{
    all_counties = in_counties;
}

std::vector<Farm*> County::get_farms(Farm_type* ft)
{
    return farms_by_type[ft];
}

std::vector<Farm_type*> County::get_farm_types_present()
{
    std::vector<Farm_type*> ftypes;
    for(auto& ft_farms_pair : farms_by_type)
    {
        ftypes.push_back(ft_farms_pair.first);
    }
    return ftypes;
}


double County::get_ocov_weight(Farm_type* ft)
{
    return ocov_weights[ft];
}

double County::get_dcov_weight(Farm_type* ft)
{
    return dcov_weights[ft];
}

double County::get_o_farm_weight_sum(Farm_type* ft)
{
    return o_farm_weight_norms.at(ft);
}

double County::get_d_farm_weight_sum(Farm_type* ft)
{
    return d_farm_weight_norms.at(ft);
}

double County::get_unnormed_o_weight(Farm_type* ft)
{
    return this->get_ocov_weight(ft) * double(this->get_n_farms(ft));
}

double County::get_normed_o_weight(Farm_type* ft)
{
    return normed_o_weight.at(ft);
}

/// Generates one shipment originating from this county.
County* County::get_shipment_destination(Farm_type* ft)
{
    County* destination = nullptr;
    if(!is_set_shipment) //If the shipping prob has not been created for this county before
    {
        this->update_shipping_probabilities(*all_counties);
    }

    if(farms_by_type.find(ft) == farms_by_type.end())
    {
        std::cout << "There are no farms of type " << ft->get_species() <<
                     " in county " << id << ". Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t ft_i = ft->get_index();
    size_t k = n_outcomes.at(ft_i);
    unsigned int outcomes[k];
    gsl_ran_multinomial(R, k, 1, shipping_probabilities.at(ft_i), outcomes);
    for(size_t i = 0; i < k; i++)
    {
        if(outcomes[i] == 1)
        {
            destination = shipping_outcomes->at(i);
            break;
        }
    }
    if(destination != nullptr)
        return destination;
    else
    {
        std::cout << "Failed to create a shipping destination." << std::endl;
        std::cout << "Exiting...";
        exit(EXIT_FAILURE);
    }
    return destination;
}

double County::cov_weight_fun(std::vector<double> cov_values,
                              std::vector<double> cov_parameters)
{
    double exponent = 0.0;
    for(size_t i = 0; i<cov_values.size(); i++)
    {
        exponent += cov_values.at(i) * cov_parameters.at(i);
    }
    return std::exp(exponent);
}


void County::calculate_centroid()
{
    if(member_farms.size() > 0)
    {
        double x_mean, y_mean;
        double x_sum = 0.0;
        double y_sum = 0.0;

        for(Farm* f : member_farms)
        {
            x_sum += f->get_x();
            y_sum += f->get_y();
        }

        x_mean = x_sum / member_farms.size();
        y_mean = y_sum / member_farms.size();
        this->set_position(x_mean, y_mean);
    }
    else
    {
        std::cout << "Error: The centroid of county " << get_id()
                  << "cannot be estimated because there are no farms in the county."
                  << " Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void County::print_bools()
{
	std::cout << "For " << id << std::endl;
	std::cout << "id = " << is_set_id << std::endl;
	std::cout << "position = " << is_set_position << std::endl;
	std::cout << "area = " << is_set_area << std::endl;
	std::cout << "state = " << is_set_state << std::endl;
	std::cout << "region = " << region_initialized << std::endl;
}

void County::set_initialized(bool& parameter)
{
    parameter = true;
    all_initialized();
}

void County::all_initialized()
{
    Region::all_initialized();
//    if(is_set_area and is_set_state and
//       is_set_id and region_initialized)
    if(is_set_area and
       is_set_id and region_initialized)
    {
        county_initialized = true;
    }
}

bool County::is_initialized()
{
	return (is_set_area and is_set_id and region_initialized);
}
