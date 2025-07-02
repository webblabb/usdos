//© 2024 Colorado State University
#include "County.h"
#include "Farm.h"
#include "State.h"
#include "Shipment_kernel.h"
#include "USAMMv2_parameters.h"
#include "USAMMv3_parameters.h"
#include "shared_functions.h"

#include <limits>
#include <algorithm>
#include <utility>
#include <gsl/gsl_randist.h>

County::County(std::string id, std::string kernel_str, int idx) :
    Region(id), kernel_str(kernel_str), idx(idx), area(0.0)
{
    verbose = verboseLevel;
    fips_code = std::stoi(id);
    type = "county";
    size_t seed = generate_distribution_seed();
    R = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(R, seed);
    tot_n_mkts = 0;
    tot_mkt_volume = 0.0;
}

County::County(std::string id, double x, double y, std::string kernel_str) :
    Region(id, x, y), kernel_str(kernel_str), area(0.0)
{
    fips_code = std::stoi(id);
    type = "county";
    tot_n_mkts = 0;
    tot_mkt_volume = 0.0;
}

County::~County()
{
    for(size_t i = 0; i < shipping_probabilities.size(); i++)
    {
        delete shipping_probabilities.at(i);
    }
    gsl_rng_free(R);

    for(auto it : receiver_prob_distributions_by_type)
    {
        gsl_ran_discrete_t* t = it.second;
        gsl_ran_discrete_free(t);
    }

    if(slaughter_facility_lookup_table != nullptr)
    {
        gsl_ran_discrete_free(slaughter_facility_lookup_table);
    }
}

//Measures the distance to all counties and calculates probabilities to send to them.
//Arguments: a vector of all th counties, a pointer to a function describing the
//kernel and a function that calculates the distance between two point objects (pointers).
void County::update_shipping_probabilities_USAMMv2(std::vector<County*>& in_counties)
{
    for(size_t i = 0; i < shipping_probabilities.size(); i++)
    {
        delete [] shipping_probabilities.at(i);
    }

    if(!is_initialized())
    {
        std::cout << "Make sure all of the following are set before attempting to "
                  << "initialize shipment probabilities." << std::endl;
        print_bools();
        not_initialized();
    }

    shipping_probabilities.resize(all_farm_types.size());
    n_outcomes.resize(all_farm_types.size());
    shipping_outcomes = &in_counties;

    //Calculate shipping probabilities for all farm-types (species).
    for(Farm_type* current_ft : all_farm_types)
    {
        size_t ft_i = current_ft->Farm_type::get_index();
        size_t n = in_counties.size();
        n_outcomes.at(ft_i) = n;
        shipping_probabilities.at(ft_i) = new double[n];

        std::vector<double> probabilities;
        probabilities.reserve(in_counties.size());
        //The kernel is based on the states current a and b and gives the shape of the distance dependence.
        Shipment_kernel kCVI(this->get_parent_state()->get_aCVI(current_ft),
                             this->get_parent_state()->get_bCVI(current_ft),
                             kernel_str, true);
        double normalization_sum = 0.0;

        //Get all kernel values and keep track of the total for use when normalizing.
        for(auto c : in_counties)
        {
            //Destination inflow weight based on number of possible farms that can receive in the receiving county,
            //therefore n-1 if sending within the county itself.
            double d_farm_weight = c->get_d_unnormalized_prem_weight_sum(current_ft);
            if(c == this)
            {
                d_farm_weight -= this->weighted_avg_d_prem_weights.at(current_ft);
            }

            //Inflow is s * sum of all farms individual weights.
            double inflow_weight = c->get_parent_state()->get_s(current_ft) * d_farm_weight;

            //Kernel value * inflow of state of destination county * sum of individual farms' weight in destination county.
            double kernel_value = kCVI.kernel(this, c);
            double unnormalized_probability = kernel_value * inflow_weight;

            if(unnormalized_probability == 0.0 and
               inflow_weight != 0.0)
            {
                unnormalized_probability = std::numeric_limits<double>::min();
            }
            if(std::isinf(unnormalized_probability) or std::isnan(unnormalized_probability))
            {
                std::cout << "Shipping probability for county " << this->id << " is "
                          << unnormalized_probability << ". Exiting..." << std::endl;
                exit(EXIT_FAILURE);
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
//void County::set_farms(const std::vector<Farm*>& in_farms)
//{
//    for(Farm* in_farm : in_farms)
//    {
//        this->add_premises(in_farm);
//    }
//}

//Adds one single farm that belongs to this county by passing a pointer to it.
//Also called internally by set_farms
void County::add_premises(Farm* in_farm, const std::vector<Farm_type*>& in_all_farm_types,
                          std::map<Farm_type*, std::set<Prem_class*>> in_all_prem_classes)
{
    all_farm_types = in_all_farm_types;
    member_premises.push_back(in_farm);
    Prem_class* in_pcl = in_farm->get_prem_class();
    std::string pcl_tag = in_pcl->tag;
    size_t pcl_idx = size_t(in_pcl->idx);

    all_premises_by_class[pcl_tag].push_back(in_farm);
    Farm_type* this_ft = in_farm->get_farm_type();
    in_farm->set_parent_county(this);
    in_farm->set_idx_in_county(int(member_premises.size()-1));

    if(n_premises_of_type_class.empty())
    {
        n_premises_of_type_class.resize(all_farm_types.size());
        for(Farm_type* ft : in_all_farm_types)
        {
            n_premises_of_type_class[ft->get_index()].resize(in_all_prem_classes.at(ft).size(), 0);
        }
    }
    n_premises_of_type_class[this_ft->get_index()][pcl_idx] += 1;

    if(pcl_tag == "Fdl")
    {
        for(Farm_type* ft : all_farm_types) //Add feedlots to all types, so they contribute to shipments for all types.
        {
            farms_by_type[ft].push_back(in_farm);
            if(n_feedl_by_type.find(ft) == n_feedl_by_type.end())
            {
                n_feedl_by_type[ft] = 0;
            }
            n_feedl_by_type.at(ft) += 1;
        }
    }
    else if(pcl_tag == "Mkt")
    {
        for(Farm_type* ft : all_farm_types)
        {
            farms_by_type[ft].push_back(in_farm); //Add markets to all types, so they contribute to shipments for all types.
        }
        this->add_market(this_ft, in_farm->get_USAMMv3_unbinned_size()*52); //Scale up weekly vol to yearly. Use quarter 3 (index 2) for USAMM related sizes.
    }
    else
    {
        farms_by_type[this_ft].push_back(in_farm);
        if(n_farms_by_type.find(this_ft) == n_farms_by_type.end())
        {
            n_farms_by_type[this_ft] = 0;
        }
        n_farms_by_type.at(this_ft) += 1;
    }
}

void County::add_market(Farm_type* ft, double vol)
{
    if(mkt_volumes.find(ft) == mkt_volumes.end())
    {
        mkt_volumes[ft] = 0.0;
    }
    mkt_volumes.at(ft) += vol;
    tot_mkt_volume += vol;
    if(n_mkt_by_type.find(ft) == n_mkt_by_type.end())
    {
        n_mkt_by_type[ft] = 0;
    }
    n_mkt_by_type.at(ft) += 1;
    ++tot_n_mkts;
}

void County::set_area(double in_area)
{
    area = in_area;
    set_initialized(is_set_area);
}

double County::pop_covariate(std::string cov_name, std::vector<std::string>& cov_names,
                             std::vector<double>& cov_values)
{
    size_t covariate_index = 0;
    bool covariate_found = false;
    for(; covariate_index < cov_names.size(); covariate_index++)
    {
        if(cov_names.at(covariate_index).compare(cov_name) == 0)
        {
            //MarketVolume covariate found.
            covariate_found = true;
            break;
        }
    }

    double covariate_value = 0.0;
    if(covariate_found)
    {
        //Save the market covariate
        covariate_value = cov_values.at(covariate_index);
        //Delete the market covariate from the vectors of other covariates.
        std::swap(cov_names.at(covariate_index), cov_names.back());
        cov_names.pop_back();
        std::swap(cov_values.at(covariate_index), cov_values.back());
        cov_values.pop_back();
    }

    return covariate_value;
}

void County::set_covariates_USAMMv2(std::map<Farm_type*, USAMMv2_parameters>& up_map)
{
    for(auto& name_up_pair : up_map)
    {
        Farm_type* up_ft = name_up_pair.second.get_farm_type();
        county_ocov_names[up_ft] = name_up_pair.second.get_county_ocov_names();
        ocov_values[up_ft] = name_up_pair.second.get_county_o_covs(this);
        county_dcov_names[up_ft] = name_up_pair.second.get_county_dcov_names();
        dcov_values[up_ft] = name_up_pair.second.get_county_d_covs(this);
    }
}

void County::set_covariates_USAMMv3(std::vector<USAMMv3_parameters>& up_vec)
{
    for(auto& up : up_vec)
    {
        Farm_type* up_ft = up.getFarmType();
        county_ocov_names[up_ft] = up.get_county_ocov_names();
        ocov_values[up_ft] = up.get_county_o_covs(this);
        county_dcov_names[up_ft] = up.get_county_dcov_names();
        dcov_values[up_ft] = up.get_county_d_covs(this);
    }
}


void County::update_covariate_weights_USAMMv2(std::map<Farm_type*, USAMMv2_parameters>& up_map,
                                              std::string time_period)
{
    county_dcov_weights.clear();
    county_dcov_weights.resize(up_map.size());
    for(auto& name_up_pair : up_map)
    {
        USAMMv2_parameters& up = name_up_pair.second;
        Farm_type* up_ft = up.get_farm_type();

        county_ocov_weights[up_ft] = 0.0;
        county_dcov_weights[up_ft->get_index()] = 0.0;
        o_unnormalized_prem_weight_sum[up_ft] = 0.0;
        d_unnormalized_prem_weight_sum[up_ft] = 0.0;
//        omarket_cov_weights[up_ft] = 0.0;
//        dmarket_cov_weights[up_ft] = 0.0;
        weighted_avg_d_prem_weights[up_ft] = 0.0;
        receiver_weights[up_ft].clear();
        receiver_weights[up_ft].reserve(farms_by_type.at(up_ft).size());
        receiver_weights_corr_farms[up_ft].clear();
        receiver_weights_corr_farms[up_ft].reserve(farms_by_type.at(up_ft).size());

        if(farms_by_type.at(up_ft).size() > 0 )
        {
            //Origin county-level covariates.
            std::vector<double> county_ocov_pars(county_ocov_names.at(up_ft).size(), 0.0);
            for(size_t i = 0; i< county_ocov_pars.size(); i++)
            {
                county_ocov_pars[i] = up.get_ocov_county_par(county_ocov_names.at(up_ft).at(i), time_period);
            }
            county_ocov_weights[up_ft] = cov_weight_fun(ocov_values.at(up_ft), county_ocov_pars);

            //Destination county-level covariates.
            std::vector<double> county_dcov_pars(county_dcov_names.at(up_ft).size(), 0.0);
            for(size_t i = 0; i< county_dcov_pars.size(); i++)
            {
                county_dcov_pars[i] = up.get_dcov_county_par(county_dcov_names.at(up_ft).at(i), time_period);
            }
            county_dcov_weights[up_ft->get_index()] = cov_weight_fun(dcov_values.at(up_ft), county_dcov_pars);

            double temp_oweight_sum = 0.0;
            double temp_dweight_sum = 0.0;

            //Exponents defaults to zero since if d/o covariates in USAMM are turned off, the o/d weights will be equal to number of premises in county (sum(1.0n^0) over all premises in county)
            double o_farm_size_cov_exp = 0.0;
            double o_farm_size_cov_coeff = 1.0;
            double d_farm_size_cov_exp = 0.0;
            double d_farm_size_cov_coeff = 1.0;
            double o_feedl_size_cov_exp = 0.0;
            double o_feedl_size_cov_coeff = 1.0;
            double d_feedl_size_cov_exp = 0.0;
            double d_feedl_size_cov_coeff = 1.0;
            double o_mkt_size_cov_exp = 0.0;
            double o_mkt_size_cov_coeff = 1.0;
            double d_mkt_size_cov_exp = 1.0;
            double d_mkt_size_cov_coeff = 1.0;

            //Get the farm-size usamm parameters.
            if(up.has_ocov_frm_par())
            {
                o_farm_size_cov_exp = up.get_ocov_prem_par("FarmExp", time_period);
            }
            if(up.has_dcov_frm_par())
            {
                d_farm_size_cov_exp = up.get_dcov_prem_par("FarmExp", time_period);
            }
//            o_farm_size_cov_coeff = std::exp(up.get_ocov_prem_par("FarmCoeff", time_period));
//            d_farm_size_cov_coeff = std::exp(up.get_dcov_prem_par("FarmCoeff", time_period));

            //Get the feedlot size parameters if they are present, otherwise it's kept at zero and feedlot weight = feedlot n animals.
            if(up.has_ocov_feedl_par())
            {
                o_feedl_size_cov_exp = up.get_ocov_prem_par("FeedlotExp", time_period);
                o_feedl_size_cov_coeff = std::exp(up.get_ocov_prem_par("FeedlotCoeff", time_period));
            }
            if(up.has_dcov_feedl_par())
            {
                d_feedl_size_cov_exp = up.get_dcov_prem_par("FeedlotExp", time_period);
                d_feedl_size_cov_coeff = std::exp(up.get_dcov_prem_par("FeedlotCoeff", time_period));
            }

            //Get the market size parameters if they are present, otherwise it's kept at zero and market weight = market n animals.
            if(up.has_ocov_mkt_par())
            {
                o_mkt_size_cov_exp = up.get_ocov_prem_par("MarketExp", time_period);
//                o_mkt_size_cov_exp = 1.0;
                o_mkt_size_cov_coeff = std::exp(up.get_ocov_prem_par("MarketCoeff", time_period));
            }
            if(up.has_dcov_mkt_par())
            {
                d_mkt_size_cov_exp = up.get_dcov_prem_par("MarketExp", time_period);
//                d_mkt_size_cov_exp = 1.0;
                d_mkt_size_cov_coeff = std::exp(up.get_dcov_prem_par("MarketCoeff", time_period));
            }

            //Calculate the premises' weight norms (sums) and set premises' weights.
            double temp_avg_prem_weight = 0.0;
            for(Farm* f : farms_by_type.at(up_ft))
            {
                double ow = 0.0;
                double dw = 0.0;
                double n_animals = 0.0;
                if(f->is_market())
                {
                    n_animals = double(f->get_size_allSpecies());
                }
                else
                {
                    n_animals = double(f->get_USAMMv3_unbinned_size());
                }
                if(n_animals > 0.0)
                {
                    if(f->is_farm())
                    {
                        ow = o_farm_size_cov_coeff * std::pow(n_animals / national_avg_farm_vol.at(up_ft),
                                                              o_farm_size_cov_exp);
                        dw = d_farm_size_cov_coeff * std::pow(n_animals / national_avg_farm_vol.at(up_ft),
                                                              d_farm_size_cov_exp);
                    }
                    else if(f->is_feedlot())
                    {
                        ow = o_feedl_size_cov_coeff * std::pow(n_animals / national_avg_feedl_vol.at(up_ft),
                                                               o_feedl_size_cov_exp);
                        dw = d_feedl_size_cov_coeff * std::pow(n_animals / national_avg_feedl_vol.at(up_ft),
                                                               d_feedl_size_cov_exp);
                    }
                    else if(f->is_market())
                    {
                        ow = o_mkt_size_cov_coeff * std::pow(n_animals / national_avg_mkt_vol.at(up_ft),
                                   o_mkt_size_cov_exp);
                        dw = d_mkt_size_cov_coeff * std::pow(n_animals / national_avg_mkt_vol.at(up_ft),
                                                               d_mkt_size_cov_exp);
                    }
                    else
                    {
                        ow = o_mkt_size_cov_coeff * std::pow(52 * n_animals / national_avg_mkt_vol.at(up_ft),
                                                             o_mkt_size_cov_exp);
                        dw = d_mkt_size_cov_coeff * std::pow(52 * n_animals / national_avg_mkt_vol.at(up_ft),
                                                             d_mkt_size_cov_exp);
                    }

                    if(std::isinf(ow) or std::isnan(ow) or
                       std::isinf(dw) or std::isnan(dw))
                    {
                        std::cout << "The premises weights of premises with id " << f->get_id() <<
                                     " in county " << this->id << " is " << ow << " (origin) and " <<
                                     dw << " (destination). Exiting." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                    ow = ow * county_ocov_weights.at(up_ft); //Weigh individual premises weight by county-level origin covariate weight.
                    dw = dw * county_dcov_weights.at(up_ft->get_index());
                }
                f->set_unnormalized_oweight(ow, up_ft);
                f->set_unnormalized_dweight(dw, up_ft);
                temp_oweight_sum += ow;
                temp_dweight_sum += dw;
                temp_avg_prem_weight += ow*dw;
                receiver_weights.at(up_ft).push_back(dw);
                receiver_weights_corr_farms.at(up_ft).push_back(f);

            }

            o_unnormalized_prem_weight_sum[up_ft] = temp_oweight_sum;
            d_unnormalized_prem_weight_sum[up_ft] = temp_dweight_sum;

            //Update the average premises size (incl markets) used for determining the probability of sending
            //a shipment within the county adjusted for the fact that the origin premises does not send to itself.
            if(temp_oweight_sum > 0.0)
            {
                weighted_avg_d_prem_weights[up_ft] = temp_avg_prem_weight / temp_oweight_sum;
            }

            //Initiate the gsl discrete probability distribution object used to
            //generate shipment receivers for this county. Delete the old one 1st
            //to avoid memory leak.
            if(receiver_prob_distributions_by_type.find(up_ft) !=
               receiver_prob_distributions_by_type.end())
            {
                gsl_ran_discrete_t* t = receiver_prob_distributions_by_type.at(up_ft);
                gsl_ran_discrete_free(t);
            }
            receiver_prob_distributions_by_type[up_ft] =
                gsl_ran_discrete_preproc(receiver_weights.at(up_ft).size(), //Number of possible outcomes
                                         receiver_weights.at(up_ft).data()); //Data returns a pointer to the internal array of probabilities.
        }
    }
}

void County::update_covariate_weights_USAMMv3(std::vector<USAMMv3_parameters>& up_vec,
                                              std::string time_period)
{
    county_dcov_weights.clear();
    county_dcov_weights.resize(up_vec.size());
    for(auto& up : up_vec)
    {
        Farm_type* up_ft = up.getFarmType();

        county_ocov_weights[up_ft] = 0.0;
        county_dcov_weights[up_ft->get_index()] = 0.0;

        //Origin county-level covariates.
        std::vector<double> county_ocov_pars(county_ocov_names.at(up_ft).size(), 0.0);
        for(size_t i = 0; i< county_ocov_pars.size(); i++)
        {
            county_ocov_pars[i] = up.get_ocov(county_ocov_names.at(up_ft).at(i), time_period);
        }
        county_ocov_weights[up_ft] = cov_weight_fun(ocov_values.at(up_ft), county_ocov_pars);

        //Destination county-level covariates.
        std::vector<double> county_dcov_pars(county_dcov_names.at(up_ft).size(), 0.0);
        for(size_t i = 0; i< county_dcov_pars.size(); i++)
        {
            county_dcov_pars[i] = up.get_dcov(county_dcov_names.at(up_ft).at(i), time_period);
        }
        county_dcov_weights[up_ft->get_index()] = cov_weight_fun(dcov_values.at(up_ft), county_dcov_pars);
    }
}

void County::init_USAMMv3_shipment_vectors(std::vector<USAMMv3_parameters>& up_vec)
{
    if(prem_classes_by_type.empty())
    {
        prem_classes_by_type.resize(up_vec.size());
        for(auto& up : up_vec)
        {
            Farm_type* up_ft = up.getFarmType();
            for(auto& pcl : up.getPremClasses())
            {
                if(n_premises_of_type_class[up_ft->get_index()][pcl->idx] > 0)
                {
                    prem_classes_by_type[up_ft->get_index()].push_back(pcl);
                }
            }
        }
    }

    //Premsize representation
//    dest_premsize_rep_CVI_by_quarter.clear();
//    origin_premsize_rep_CVI_by_quarter.clear();
    n_time_periods = up_vec.begin()->get_n_time_periods();
//    size_t n_ft = up_vec.size();
//    dest_premsize_rep_CVI_by_quarter.resize(n_time_periods, Vec_d_4d(n_ft));
//    origin_premsize_rep_CVI_by_quarter.resize(n_time_periods, Vec_d_4d(n_ft));

    //N-vector
    if(N_vec_CVI.empty())
    {
        N_vec_CVI.resize(up_vec.size());
    }

    for(USAMMv3_parameters& up : up_vec)
    {
        Farm_type* ft = up.getFarmType();


        Vec_i_2d& N_vec_CVI_ft = N_vec_CVI.at(ft->get_index());
        if(up.getUSAMMv3SpeciesCode() != "s" and !N_vec_CVI_ft.empty())
        {
            //When simulating cattle, the N_vec doesn't change with new USAMM parameters and
            //since it's already been initialized we can skip this step. If simulating swine,
            //the N_vec changes based on which prems are associated with COM ids each rep.
            continue;
        }

        std::vector<Prem_class*> prem_classes = up.getPremClasses();
        N_vec_CVI_ft.clear();
        N_vec_CVI_ft.resize(prem_classes.size());
        for(Prem_class* pcl : prem_classes)
        {
            std::set<int> size_bin_set = up.getSizeBins(ft, pcl);
            std::vector<int> bins(size_bin_set.begin(), size_bin_set.end());
            std::vector<int>& N_vec_CVI_ft_pcl = N_vec_CVI_ft.at(pcl->idx);
            N_vec_CVI_ft_pcl.resize(bins.size(), 0);

            std::map<int, int> binsize_to_index;
            for(size_t i=0; i<bins.size(); ++i)
            {
                binsize_to_index[bins[i]] = i;
            }

            if(pcl->tag == "Fdl" or pcl->tag == "Mkt")
            {
                for(Farm* f : member_premises)
                {
                    if(f->get_prem_class() == pcl and f->getCOMIdentifier() == nullptr)
                    {
                        size_t bin_idx = binsize_to_index.at(f->get_USAMMv3_binned_size());
                        ++N_vec_CVI_ft[pcl->idx][bin_idx];
                    }
                }
            }
            else
            {
                for(Farm* f : member_premises)
                {
                    if(f->get_prem_class() == pcl and
                       f->get_farm_type() == ft and
                       f->getCOMIdentifier() == nullptr)
                    {
                        size_t bin_idx = binsize_to_index.at(f->get_USAMMv3_binned_size());
                        ++N_vec_CVI_ft[pcl->idx][bin_idx];
                    }
                }
            }
        }
    }
}

void County::init_ft_pcl_vec(const std::map<Farm_type*, std::set<Prem_class*>>& ft_pcl_map)
{
    all_farm_types.resize(ft_pcl_map.size());
    d_CVI_lambda_sums.resize(ft_pcl_map.size());
    for(const auto ft_pcl_pair : ft_pcl_map)
    {
        Farm_type* ft = ft_pcl_pair.first;
        all_farm_types[ft->get_index()] = ft;
        d_CVI_lambda_sums[ft->get_index()] = std::vector<double>(ft_pcl_pair.second.size(), -1.0);
    }
}

//void County::update_cc_CVI_shipment_rate(std::vector<USAMMv3_parameters>& up_vec,
//                                         Farm_type* ft, std::string time_period)
//{
//    //If it's the first time here, or shipment rates have been reset due to new USAMMv3 parameters,
//    //are used (new time period). Reshape the premsize rep to the correct dimensions.
//    if(cc_CVI_shipment_rates_by_quarter.find(time_period) == cc_CVI_shipment_rates_by_quarter.end())
//    {
//        cc_CVI_shipment_rates_by_quarter[time_period] = new Vec_d_4d;
//        c_origin_CVI_shipment_rates_by_quarter[time_period] = new Vec_d_3d;
//    }
//    Vec_d_4d& cc_shipment_rates = (*cc_CVI_shipment_rates_by_quarter[time_period]);
//    Vec_d_3d& c_origin_shipment_rates = (*c_origin_CVI_shipment_rates_by_quarter[time_period]);
//
//    if(cc_shipment_rates.empty()) //If this is empty, c_origin_shipment_rates will also be empty as they are always build together.
//    {
//        cc_shipment_rates.resize(up_vec.size());
//        c_origin_shipment_rates.resize(up_vec.size());
//    }
//
//    //For each farm type (beef/dairy) create one set of premsize rep for every possible
//    //sending prem class (farm, feedlot, market).
//    int up_ft_idx = ft->get_index();
//    if(cc_shipment_rates[up_ft_idx].empty())
//    {
//        USAMMv3_parameters& up = up_vec[ft->get_index()];
//        Vec_d_3d& cc_rates_ft = cc_shipment_rates[up_ft_idx];
//        Vec_d_2d& c_origin_rates_ft = c_origin_shipment_rates[up_ft_idx];
//        std::vector<Prem_class*> up_prem_classes = up.getPremClasses();
//        int n_classes = int(up_prem_classes.size());
//
//        cc_rates_ft.resize(n_classes, Vec_d_2d(n_classes,
//                                               std::vector<double>(all_counties->size(), 0.0)));
//        c_origin_rates_ft.resize(n_classes, Vec_d_1d(n_classes, 0.0));
//
//        //Prepare some variables to reduce the scope and improve calc speed.
//        double outflow = up.get_outflowCVI(this->parent_state->get_code(), time_period);
//        double ocov_weight = this->get_county_ocov_weight(ft);
//
//        //At this point cc_rates_ft, which is a reference to the premsize rep for a given Farm_type,
//        //is fully initialized with indices: sending prem class idx, receiving prem class idx,
//        //receiving county idx. The last index matches the vector all_counties, which is a vector of
//        //pointers to counties and is shared between all counties. The values of cc_shipment_rates is the
//        //summed total prem-to-prem shipment rate for all premises of the relevant types in the two
//        //counties (this county sending to the county indicated by the last idex).
//        //In order to reduce calculation in the shipment generation process the inflow and outflow
//        //weights of respective counties state as well as their origin and destination covariate
//        //weights are included here in addition to the premises sizes.
//        for(Prem_class* o_pcl : up.getPremClasses())
//        {
//            //These are the bin sizes v evaluated as v^phi_o
//            std::vector<double> o_bin_weights = up.get_evaluated_origin_bin_weights(o_pcl, time_period);
//            std::vector<int>& N_vec_CVI_ft_pcl = N_vec_CVI.at(ft->get_index()).at(o_pcl->idx); //This is the number of premises of this class in this county.
//            double o_lambda_sum = 0.0;
//            for(size_t i=0; i<o_bin_weights.size(); ++i)
//            {
//                o_lambda_sum += N_vec_CVI_ft_pcl[i] * o_bin_weights[i]; //Number of prems of spec. bin * their individual origin rate or weight.
//            }
//
//            //Calculate the total lambda of all prems of this prem class sending to every other
//            //combination of prem class and county.
//            for(Prem_class* d_pcl : up.getPremClasses())
//            {
//                double pcl_pcl_c = up.get_c(o_pcl, d_pcl, time_period);
//
//                if(o_lambda_sum == 0)
//                {
//                    cc_rates_ft[o_pcl->idx][d_pcl->idx] = std::vector<double>(all_counties->size(), 0.0);
//                }
//                else
//                {
//                    //I put the kernel here mostly since we might end up with prem-class specific
//                    //kernel parameters at some point.
//                    Shipment_kernel kCVI(this->get_parent_state()->get_aCVI(ft),
//                                         this->get_parent_state()->get_bCVI(ft),
//                                         kernel_str, true);
//
//                    for(size_t c_idx=0; c_idx<all_counties->size(); ++c_idx)
//                    {
//                        County* c = all_counties->at(c_idx);
//                        double d_lambda_sum = c->get_d_CVI_shipment_rate_sum(d_pcl, time_period, up);
//                        double cc_lambda = 0.0;
//                        if(d_lambda_sum > 0.0)
//                        {
//                            double k_val_CVI = kCVI.kernel(this, c);
//                            cc_lambda = o_lambda_sum *
//                                        outflow *
//                                        ocov_weight *
//                                        d_lambda_sum *
//                                        up.get_inflowCVI(c->parent_state->get_code(), time_period) *
//                                        c->get_county_dcov_weight(ft) *
//                                        pcl_pcl_c *
//                                        k_val_CVI;
//                        }
//                        cc_rates_ft[o_pcl->idx][d_pcl->idx][c_idx] = cc_lambda;
//                        c_origin_rates_ft[o_pcl->idx][d_pcl->idx] += cc_lambda;
//                    }
//                }
//            }
//        }
//    }
//}

void County::unset_shipping_probabilities(bool full_reset)
{
    size_t n_ft = all_farm_types.size();
    size_t n_prem_classes = d_CVI_lambda_sums[0].size();

    if(full_reset)
    {
//        dest_premsize_rep_CVI_by_quarter.clear();
//        dest_premsize_rep_CVI_by_quarter.resize(n_time_periods, Vec_d_4d(n_ft));
//        origin_premsize_rep_CVI_by_quarter.clear();
//        origin_premsize_rep_CVI_by_quarter.resize(n_time_periods, Vec_d_4d(n_ft));
        weighted_o_N_vec_CVI.clear();
        weighted_o_N_vec_CVI.resize(n_time_periods, Vec_d_3d(n_ft));
        weighted_d_N_vec_CVI.clear();
        weighted_d_N_vec_CVI.resize(n_time_periods, Vec_d_3d(n_ft));
        within_county_adj_terms_CVI.clear();
        within_county_adj_terms_CVI.resize(n_time_periods, Vec_d_2d(n_ft));
        outgoing_county_CVI_weight.clear();
        outgoing_county_CVI_weight.resize(n_ft, Vec_d_2d(n_prem_classes, Vec_d_1d(n_prem_classes, 0.0)));

//        for(auto tp_vec_pair : cc_CVI_shipment_rates_by_quarter)
//        {
//            delete tp_vec_pair.second;
//        }
//        cc_CVI_shipment_rates_by_quarter.clear();

//        for(auto tp_vec_pair : c_origin_CVI_shipment_rates_by_quarter)
//        {
//            delete tp_vec_pair.second;
//        }
//        c_origin_CVI_shipment_rates_by_quarter.clear();
    }

    is_set_shipment = false;
    internal_oprem_weights_CVI.clear();
    internal_oprem_weight_sums_CVI.clear();
    internal_dprem_weights_CVI.clear();
    d_CVI_lambda_sums.clear();
    d_CVI_lambda_sums.resize(all_farm_types.size(),
                         std::vector<double>(n_prem_classes, -1.0));
}

void County::normalize_shipping_weights(Farm_type* ft, double o_norm)
{
    for(Farm* f : farms_by_type.at(ft))
    {
        f->set_normalized_oweight(f->get_unnormalized_oweight(ft) / o_norm, ft);
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

void County::set_national_avg_farm_vol(Farm_type* ft, double vol)
{
    national_avg_farm_vol[ft] = vol;
}

void County::set_national_avg_feedl_vol(Farm_type* ft, double vol)
{
    national_avg_feedl_vol[ft] = vol;
}

void County::set_national_avg_mkt_vol(Farm_type* ft, double vol)
{
    national_avg_mkt_vol[ft] = vol;
}

void County::set_outgoing_county_CVI_weight(size_t fty_idx, size_t o_pcl_idx,
                                            size_t d_pcl_idx, double v)
{
    outgoing_county_CVI_weight.at(fty_idx).at(o_pcl_idx).at(d_pcl_idx) = v;
}

int County::get_n_farms(Farm_type* ft)
{
    if(n_farms_by_type.find(ft) != n_farms_by_type.end())
    {
        return n_farms_by_type.at(ft);
    }
    return 0;
}

int County::get_n_feedl(Farm_type* ft)
{
    if(n_feedl_by_type.find(ft) != n_feedl_by_type.end())
    {
        return n_feedl_by_type.at(ft);
    }
    return 0;
}

int County::get_n_mkt()
{
    return tot_n_mkts;
}

double County::get_total_farm_vol(Farm_type* ft)
{
    double sum = 0.0;
    for(Farm* f : this->get_premises(ft))
    {
        if(!f->is_feedlot() && !f->is_market())
        {
            sum += f->get_size_allSpecies();
        }
    }
    return sum;
}

double County::get_total_feedl_vol(Farm_type* ft)
{
    double sum = 0.0;
    for(Farm* f : this->get_premises(ft))
    {
        if(f->is_feedlot())
        {
            sum += f->get_size_allSpecies();
        }
    }
    return sum;
}

double County::get_total_mkt_vol()
{
    return tot_mkt_volume;
}

std::vector<Farm*>& County::get_premises(Farm_type* ft)
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

bool County::has_premises_of_type_class(Farm_type* ft, Prem_class* pcl)
{
    return n_premises_of_type_class[ft->get_index()][pcl->idx] > 0;
}

double County::get_county_ocov_weight(Farm_type* ft)
{
    return county_ocov_weights[ft];
}

double County::get_county_dcov_weight(Farm_type* ft)
{
    return county_dcov_weights[ft->get_index()];
}

double County::get_o_unnormalized_prem_weight_sum(Farm_type* ft)
{
    return o_unnormalized_prem_weight_sum.at(ft);
}

double County::get_d_unnormalized_prem_weight_sum(Farm_type* ft)
{
    return d_unnormalized_prem_weight_sum.at(ft);
}

double County::get_outgoing_county_CVI_weight(size_t fty_idx, size_t o_pcl_idx,
                                              size_t d_pcl_idx)
{
    return outgoing_county_CVI_weight.at(fty_idx).at(o_pcl_idx).at(d_pcl_idx);
}

void County::set_slaughter_probs(std::vector<double> probabilties, const std::vector<int>& facility_ids)
{
    slaughter_facility_ids = facility_ids;
    slaughter_facility_lookup_table = gsl_ran_discrete_preproc(probabilties.size(), probabilties.data());
    slaughter_probs_set = true;
}

int County::get_slaughter_destination()
{
    if(slaughter_facility_lookup_table != nullptr)
    {
        int facility_idx = int(gsl_ran_discrete(R, slaughter_facility_lookup_table));
        return slaughter_facility_ids[facility_idx];
    }
    else
    {
        return -1;
    }
}

//double County::get_d_CVI_shipment_rate_sum(Prem_class* d_pcl, std::string time_period,
//                                       USAMMv3_parameters& up)
//{
//    Farm_type* up_ft = up.getFarmType();
//    double& d_lambda_sum = d_CVI_lambda_sums[up_ft->get_index()][d_pcl->idx];
//    if(d_lambda_sum < 0.0)
//    {
//        //These are the bin sizes v evaluated as v^phi_d
//        std::vector<double> d_bin_weights = up.get_evaluated_destination_bin_weights(d_pcl, time_period);
//        const std::vector<int>& N_vec_CVI_ft_pcl = N_vec_CVI.at(up_ft->get_index()).at(d_pcl->idx); //This is the number of premises of this class in this county.
//        d_lambda_sum = 0.0;
//        for(size_t i=0; i<d_bin_weights.size(); ++i)
//        {
//            d_lambda_sum += N_vec_CVI_ft_pcl[i] * d_bin_weights[i]; //Number of prems of spec. bin * their individual origin rate or weight.
//        }
//    }
//    return d_lambda_sum;
//}

//double County::get_o_CVI_shipment_rate(Farm_type* ft, Prem_class* o_pcl,
//                                   Prem_class* d_pcl, std::string time_period,
//                                   std::vector<USAMMv3_parameters>& up_vec)
//{
//    try
//    {
//        return (*c_origin_CVI_shipment_rates_by_quarter.at(time_period)).at(ft->get_index()).at(o_pcl->idx).at(d_pcl->idx);
//    }
//    catch(const std::out_of_range& oor)
//    {
//        update_cc_CVI_shipment_rate(up_vec, ft, time_period);
//    }
//    return (*c_origin_CVI_shipment_rates_by_quarter.at(time_period))[ft->get_index()][o_pcl->idx][d_pcl->idx];
//}
//
//const std::vector<double>& County::get_d_CVI_shipment_rate_vec(Farm_type* ft, Prem_class* o_pcl,
//                                                           Prem_class* d_pcl, std::string time_period,
//                                                           std::vector<USAMMv3_parameters>& up_vec)
//{
//    try
//    {
//        return (*cc_CVI_shipment_rates_by_quarter.at(time_period)).at(ft->get_index()).at(o_pcl->idx).at(d_pcl->idx);
//    }
//    catch(const std::out_of_range& oor)
//    {
//        update_cc_CVI_shipment_rate(up_vec, ft, time_period);
//    }
//    return (*cc_CVI_shipment_rates_by_quarter.at(time_period))[ft->get_index()][o_pcl->idx][d_pcl->idx];
//}

const std::vector<double>& County::get_eval_N_vec_origin_CVI(int fty_idx, int pcl_idx,
                                                             USAMMv3_parameters& up,
                                                             int time_period_idx)
{
    Vec_d_2d& ft_weighted_o_N_vec = weighted_o_N_vec_CVI[time_period_idx][fty_idx];
    if(ft_weighted_o_N_vec.empty())
    {
        Vec_d_2d& ft_weighted_d_N_vec = weighted_d_N_vec_CVI[time_period_idx][fty_idx];
        Vec_d_1d& ft_within_county_adj_terms = within_county_adj_terms_CVI[time_period_idx][fty_idx];
        eval_N_vecs_for_fty(fty_idx, time_period_idx, up, ft_weighted_o_N_vec,
                            ft_weighted_d_N_vec, ft_within_county_adj_terms);
    }
    return ft_weighted_o_N_vec[pcl_idx];
}

const std::vector<double>& County::get_eval_N_vec_dest_CVI(int fty_idx, int pcl_idx,
                                                           USAMMv3_parameters& up,
                                                           int time_period_idx)
{
    Vec_d_2d& ft_weighted_d_N_vec = weighted_d_N_vec_CVI[time_period_idx][fty_idx];
    if(ft_weighted_d_N_vec.empty())
    {
        Vec_d_2d& ft_weighted_o_N_vec = weighted_o_N_vec_CVI[time_period_idx][fty_idx];
        Vec_d_1d& ft_within_county_adj_terms = within_county_adj_terms_CVI[time_period_idx][fty_idx];
        eval_N_vecs_for_fty(fty_idx, time_period_idx, up, ft_weighted_o_N_vec,
                            ft_weighted_d_N_vec, ft_within_county_adj_terms);
    }
    return ft_weighted_d_N_vec[pcl_idx];
}

double County::get_within_county_adj_terms_CVI(int fty_idx, int pcl_idx, USAMMv3_parameters& up,
                                               int time_period_idx)
{
    Vec_d_1d& ft_within_county_adj_terms = within_county_adj_terms_CVI[time_period_idx][fty_idx];
    if(ft_within_county_adj_terms.empty())
    {
        Vec_d_2d& ft_weighted_o_N_vec = weighted_o_N_vec_CVI[time_period_idx][fty_idx];
        Vec_d_2d& ft_weighted_d_N_vec = weighted_d_N_vec_CVI[time_period_idx][fty_idx];
        eval_N_vecs_for_fty(fty_idx, time_period_idx, up, ft_weighted_o_N_vec,
                            ft_weighted_d_N_vec, ft_within_county_adj_terms);
    }
    return ft_within_county_adj_terms[pcl_idx];
}

void County::eval_N_vecs_for_fty(int fty_idx, int time_period_idx, USAMMv3_parameters& up,
                                 Vec_d_2d& ft_weighted_o_N_vec, Vec_d_2d& ft_weighted_d_N_vec,
                                 Vec_d_1d& ft_within_county_adj_terms)
{
    const std::vector<Prem_class*>& pcl_vec = up.getPremClasses();
    size_t n_classes = pcl_vec.size();
    ft_weighted_o_N_vec.clear();
    ft_weighted_o_N_vec.resize(n_classes);
    ft_weighted_d_N_vec.clear();
    ft_weighted_d_N_vec.resize(n_classes);
    ft_within_county_adj_terms.clear();
    ft_within_county_adj_terms.resize(n_classes);
    for(Prem_class* pcl : pcl_vec)
    {
        const std::vector<double> o_bin_weights = up.get_evaluated_origin_bin_weights(pcl, time_period_idx);
        const std::vector<double> d_bin_weights = up.get_evaluated_destination_bin_weights(pcl, time_period_idx);
        const std::vector<int> N_vec = N_vec_CVI.at(fty_idx).at(pcl->idx);
        std::vector<double> o_weights(o_bin_weights.size(), 0.0);
        std::vector<double> d_weights(d_bin_weights.size(), 0.0);
        double adj_sum = 0.0;
        for(size_t i=0; i<d_bin_weights.size(); ++i)
        {
            o_weights[i] = N_vec[i] * o_bin_weights[i]; //Number of prems of spec. bin * their individual origin rate or weight.
            d_weights[i] = N_vec[i] * d_bin_weights[i]; //Number of prems of spec. bin * their individual dest rate or weight.
            adj_sum += o_bin_weights[i] * d_bin_weights[i] * N_vec[i];
        }
        ft_weighted_o_N_vec[pcl->idx].swap(o_weights);
        ft_weighted_d_N_vec[pcl->idx].swap(d_weights);
        ft_within_county_adj_terms[pcl->idx] = adj_sum;
    }
}

const std::vector<double>& County::get_USAMMv3_origin_prem_weights_CVI(size_t fty_idx, Prem_class* o_pcl,
                                                                   std::string time_period, USAMMv3_parameters& up)
{
    if(internal_oprem_weights_CVI.empty())
    {
        internal_oprem_weights_CVI.resize(all_farm_types.size(), Vec_d_2d(3)); //Max number of prem classes hard-coded to 3, change here if mor are needed.
    }

    Vec_d_1d& oweights = internal_oprem_weights_CVI[fty_idx][o_pcl->idx];
    if(oweights.empty())
    {
        oweights.resize(member_premises.size(), 0.0);
        std::map<Prem_class*, std::map<int, double>>& evaluated_bin_oweights = up.get_origin_binweight_lookup(time_period);
        if(o_pcl->tag == "Fdl" or o_pcl->tag == "Mkt")
        {
            //If fdl or market all premises of that prem class are added regardless of farm type since these classes are relevant for both beef and dairy.
            for(size_t p_idx=0; p_idx<member_premises.size(); ++p_idx)
            {
                Farm* prem = member_premises[p_idx];
                if(prem->get_prem_class() == o_pcl and prem->getCOMIdentifier() == nullptr) //Exclude COM premises.
                {
                    oweights[p_idx] = evaluated_bin_oweights.at(o_pcl).at(prem->get_USAMMv3_binned_size());
                } //Else weight stays 0.
            }
        }
        else
        {
            //The prem-class is regular farm.
            for(size_t p_idx=0; p_idx<member_premises.size(); ++p_idx)
            {
                Farm* prem = member_premises[p_idx];
                if(prem->get_prem_class() == o_pcl and prem->get_farm_type()->get_index() == fty_idx and
                   prem->getCOMIdentifier() == nullptr) //Exclude COM premises.
                {
                    oweights[p_idx] = evaluated_bin_oweights.at(o_pcl).at(prem->get_USAMMv3_binned_size());
                } //Else weight stays 0.
            }
        }
    }

    return oweights;
}

double County::get_USAMMv3_origin_prem_weight_sum_CVI(size_t fty_idx, Prem_class* o_pcl,
                                                  std::string time_period, USAMMv3_parameters& up)
{
    if(internal_oprem_weight_sums_CVI.empty())
    {
        internal_oprem_weight_sums_CVI.resize(all_farm_types.size(), Vec_d_1d(3, -1.0)); //Max number of prem classes hard-coded to 3, change here if mor are needed.
    }

    double& osum = internal_oprem_weight_sums_CVI[fty_idx][o_pcl->idx];
    if(osum < 0.0)
    {
        const Vec_d_1d& o_prem_weights = this->get_USAMMv3_origin_prem_weights_CVI(fty_idx, o_pcl, time_period, up);
        for(const double& v : o_prem_weights)
        {
            osum += v;
        }
    }
    return osum;
}

const std::vector<double>& County::get_USAMMv3_dest_prem_weights_CVI(size_t fty_idx, Prem_class* d_pcl,
                                                                 std::string time_period, USAMMv3_parameters& up)
{
    if(internal_dprem_weights_CVI.empty())
    {
        internal_dprem_weights_CVI.resize(all_farm_types.size(), Vec_d_2d(3)); //Max number of prem classes hard-coded to 3, change here if mor are needed.
    }

    Vec_d_1d& dweights = internal_dprem_weights_CVI[fty_idx][d_pcl->idx];
    if(dweights.empty())
    {
        dweights.resize(member_premises.size(), 0.0);
        std::map<Prem_class*, std::map<int, double>>& evaluated_bin_dweights = up.get_dest_binweight_lookup(time_period);
        if(d_pcl->tag == "Fdl" or d_pcl->tag == "Mkt")
        {
            //If fdl or market all premises of that prem class are added regardless of farm type since these classes are relevant for both beef and dairy.
            for(size_t p_idx=0; p_idx<member_premises.size(); ++p_idx)
            {
                Farm* prem = member_premises[p_idx];
                if(prem->get_prem_class() == d_pcl)
                {
                    dweights[p_idx] = evaluated_bin_dweights.at(d_pcl).at(prem->get_USAMMv3_binned_size());
                } //Else weight stays 0.
            }
        }
        else
        {
            //The prem-class is regular farm.
            for(size_t p_idx=0; p_idx<member_premises.size(); ++p_idx)
            {
                Farm* prem = member_premises[p_idx];
                if(prem->get_prem_class() == d_pcl and prem->get_farm_type()->get_index() == fty_idx)
                {
                    dweights[p_idx] = evaluated_bin_dweights.at(d_pcl).at(prem->get_USAMMv3_binned_size());
                } //Else weight stays 0.
            }
        }
    }
    return dweights;
}

/// Generates one shipment originating from this county.
County* County::get_shipment_destination_county(Farm_type* ft)
{
    County* destination = nullptr;
    if(!is_set_shipment) //If the shipping prob has not been created for this county before
    {
        this->update_shipping_probabilities_USAMMv2(*all_counties);
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

Farm* County::get_shipment_destination_premises(Farm_type* ft)
{
    /*
    Given that this is the destination county, selects a premises based on
    the different premises' (farms and feedlots) and market weights. If the
    destination is determined to be a market, a nullptr is returned. This works
    since at the moment multiple markets are not separated within a county but
    treated as a single total county-level market weight.
    */
    size_t outcome_idx = gsl_ran_discrete(R, receiver_prob_distributions_by_type.at(ft));
    return receiver_weights_corr_farms.at(ft).at(outcome_idx);
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
    if(member_premises.size() > 0)
    {
        double x_mean, y_mean;
        double x_sum = 0.0;
        double y_sum = 0.0;

        for(Farm* f : member_premises)
        {
            x_sum += f->get_x();
            y_sum += f->get_y();
        }

        x_mean = x_sum / member_premises.size();
        y_mean = y_sum / member_premises.size();
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

