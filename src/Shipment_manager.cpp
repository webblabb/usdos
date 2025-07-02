//© 2024 Colorado State University
#include <cmath> // floor
#include <iostream> // for error-checking output
#include <algorithm>
#include "Shipment_manager.h" // includes Farm, shared_functions
#include "Shipment_kernel.h"
#include "County.h"
#include "Status_manager.h"
#include "shared_functions.h"
#include "gsl/gsl_randist.h"

///	\param[in]	in_FIPS_map		A map of FIPS codes to farms
///	\param[in]	fipsSpMap		Sorted populations of species on farms
///	\param[in]	p	Pointer to the shared Parameters instance.
Shipment_manager::Shipment_manager(
	const std::vector<County*> in_FIPS_vec,
	const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpMap,
	const Parameters* p) :
        allCounties(in_FIPS_vec),
        fipsSpeciesMap(fipsSpMap),
        parameters(p),
        species(p->species)
{
	verbose = verboseLevel;
	initialize();
}

Shipment_manager::~Shipment_manager()
{
	gsl_rng_free(R);
}

void Shipment_manager::initialize()
{
    //Initialize the random number generator.
    size_t seed = generate_distribution_seed();
    R = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(R, seed);

    // Determine if shipping is turned off in parameters
    shipments_on = parameters->shipments_on;
    if(shipments_on){
        // copy species
        for (County* c : allCounties){ // first is FIPS id, second County*
            if(allStates_set.find(c->get_parent_state()) == allStates_set.end())
            {
                allStates_set.insert(c->get_parent_state());
            }
            for (auto& s:species){// can remove when USAMM is implemented
                if (fipsSpeciesMap->at(c->get_id()).count(s)==1){speciesFIPS[s].emplace_back(c->get_id());}
            }
        }

        //COMMUTER PROBABILITY TO OBSERVE SHIPMENT BTW STATE PAIR HARDCODED HERE.
        std::map<size_t, double> o_prop_obs_swine = { {19, 1.0},
                                                   {6, 1.0},
                                                   {26, 0.0},
                                                   {36, 1.0},
                                                   {37, 0.997028},
                                                   {55, 1.0} };
        std::map<size_t, double> d_prop_obs_swine = { {19, 0.965191},
                                                   {6, 0.0},
                                                   {26, 0.0},
                                                   {36, 0.0},
                                                   {37, 0.0},
                                                   {55, 0.0} };
        size_t max_state_idx = 56 + 1;
        p_obs_COM_lookup_table.resize(max_state_idx+1, Vec_d_1d(max_state_idx+1, 0.0));
        for(size_t ostate_id=1; ostate_id<max_state_idx+1; ++ostate_id)
        {
            double o_val = 0.0;
            if(o_prop_obs_swine.find(ostate_id) != o_prop_obs_swine.end())
            {
                o_val = o_prop_obs_swine.at(ostate_id);
            }
            for(size_t dstate_id=1; dstate_id<max_state_idx+1; ++dstate_id) //For all combinations of all states (and some that done exist but that doesn't matter).
            {
                if(ostate_id != dstate_id)
                {
                    double d_val = 0.0;
                    if(d_prop_obs_swine.find(dstate_id) != d_prop_obs_swine.end())
                    {
                        d_val = d_prop_obs_swine.at(dstate_id);
                    }
                    double v = 1.0 - (1.0 - o_val) * (1.0 - d_val);
                    p_obs_COM_lookup_table.at(ostate_id).at(dstate_id) = v;
                }
            }
        }

        if(verbose>1){std::cout << "Shipment manager constructed: "<<allCounties.size()<<" counties with premises."
                                << " Shipment kernel in use: " << parameters->shipment_kernel << "." << std::endl;}
    } else {
    	if(verbose>0){std::cout << "Shipment manager constructed but not initialized since shipments are turned off." << std::endl;}
    }
}

void Shipment_manager::makeShipmentsUSAMMv2(size_t timestep, size_t days_rem,
                                            std::string time_period, std::vector<Shipment*>& output,
                                            std::vector<Farm*>& infFarms, std::vector<Farm_type*> ft_vec)
{
    //Select what farms will be involved in the generation of shipments.
    //Running with an empty infFarms is a signal to generating a full network of shipments,
    //so then we use all farms.
    //This function does have anything to do with btb, as btb is always simulated with USAMMv3
    std::vector<Farm*> affected_farms;
    affected_farms.reserve(1000000);
    size_t day_of_year = get_day_of_year(timestep, parameters->start_day);
    if(infFarms.empty())
    {
        for(County* c : allCounties)
        {
            std::vector<Farm*> c_farms = c->get_premises();
            affected_farms.insert(affected_farms.end(), c_farms.begin(), c_farms.end());
        }
    }
    else
    {
        affected_farms.assign(infFarms.begin(), infFarms.end());
    }

    //Sort all affected farms according to farm type and state.
    std::map<Farm_type*, std::map<State*, std::vector<Farm*>>> affected_farms_by_ft_state;
    for(Farm* f : affected_farms)
    {
        if(f->is_market())
        {
            //Markets contribute to both beef and dairy shipments.
            for(Farm_type* ft : ft_vec)
            {
                affected_farms_by_ft_state[ft][f->get_parent_state()].push_back(f);
            }
        }
        else
        {
            affected_farms_by_ft_state[f->get_farm_type()][f->get_parent_state()].push_back(f);
        }
    }

    //For each farm type and state, draw the number of shipments originating from there
    //this time step from a Poisson distribution using the sum of each individual
    //farms shipping rate as the state-level rate. Then, assign these shipments to
    //the farms within the state using a multinomial distribution where each farm
    //is weighted by some parameter (here all farms have equal weight).
    for(auto& ft_and_state_farms_pair : affected_farms_by_ft_state)
    {
        Farm_type* ft = ft_and_state_farms_pair.first;
        for(auto& state_farms_pair : ft_and_state_farms_pair.second)
        {
            State* s = state_farms_pair.first;
            //Generate the number of shipments that originate from this state. This is done internally
            //within the state object based on the usamm version selected. The shipping rate within the
            //state object is expressed per day, so sum up the number of shipments occuring during the
            //timestep according to the number of days/timestep.
            int n_shipments = 0;
            for(int day_idx=0; day_idx<parameters->days_per_timestep; ++day_idx)
            {
                n_shipments += s->generate_daily_shipments(ft, days_rem);
            }

            if(n_shipments > 0)
            {
                //Construct an array to store the weights.
                size_t n_affected_farms = state_farms_pair.second.size();
                std::vector<County*> origin_counties = state_farms_pair.first->get_member_counties();
                size_t n_weight_elements = n_affected_farms + 1;
                double f_weights[n_weight_elements]; //Weight of non-infected making a shipment as last element.

                //loop over the affected farms and insert their weights into the weight array. Sum them at the same time.
                double affected_weight_sum = 0.0;
                for(size_t i = 0; i < n_affected_farms; i++)
                {
                    Farm* affected_farm = state_farms_pair.second.at(i);
                    double this_affected_weight = affected_farm->get_normalized_oweight(ft);
                    f_weights[i] = this_affected_weight;
                    affected_weight_sum += this_affected_weight;
                }

                //insert that 'other' weight after the individual affected premises' weights.
//                f_weights[n_affected_farms] = s->get_total_farm_weight(ft) - affected_weight_sum; //Last weight is the sum of the weights of all farms that are not infectious.
                f_weights[n_weight_elements-1] = 1.0 - affected_weight_sum; //Last weight is the sum of the weights of all farms that are not infectious...
                if(f_weights[n_weight_elements-1] < 0.0)
                {
                    if(f_weights[n_weight_elements-1] < -1.0e-12)
                    {
//                        std::cout << "Ship manager asdasd" << std::endl;
                    }
                    f_weights[n_weight_elements-1] = 0.0; //This is likely due to numerical problems with lots of prems with small weights. Force zero, so no negative value goes into the multinomial dist.
                }

                //Create an array to save the outcomes (number of shipment from each corresponding premises/market in the weight vector).
                unsigned int f_outcome[n_weight_elements];
                std::fill_n(f_outcome, n_weight_elements, 0); //Initiate all outcomes to 0.
                gsl_ran_multinomial(R, n_weight_elements, n_shipments, f_weights, f_outcome); //The weights are normalized internally in gsl_ran_multinomial.

                //Generate shipments from the infected farms (0 - n_affected_farms).
                for(size_t i = 0; i < n_affected_farms; i++)
                {
                    if(f_outcome[i] > 0)
                    {
                        Farm* current_farm = state_farms_pair.second[i];
                        for(size_t j = 0; j < f_outcome[i]; j++)
                        {
                            Shipment* s = generateShipmentUSAMMv2(current_farm, timestep, day_of_year, time_period);
                            output.push_back(s);
                        }
                    }
                }
            }
        }
    }
}
// 7th Sep. 2022: I checked the shipment bit by comparing networks simulated by USDOS
// against networks simulated with USAMM using the same posterior sample. Both swine
// and cattle (incl. feedlots in dairy) look fine.
void Shipment_manager::makeShipmentsUSAMMv3(size_t timestep, size_t day_of_year, std::string time_period,
                                            int time_period_idx, std::vector<Shipment*>& output,
                                            std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec)
{
    int n_ships_generated = 0;
    const Vec_d_2d& county_distance_matrix = up_vec[0].get_county_distance_matrix();
//    int delta_t = parameters->USAMM_temporal_n_timesteps[time_period_idx];
    int delta_t = 1;

    //Determine the index of the market prem class (just for convenience/speed).
    Prem_class* mkt_pcl = nullptr;
    Prem_class* fdl_pcl = nullptr;
    std::vector<Prem_class*> pcl_vec = up_vec.begin()->getPremClasses();
    for(Prem_class* pcl : pcl_vec)
    {
        if(pcl->tag == "Mkt")
        {
            mkt_pcl = pcl;
        }
        if(pcl->tag == "Fdl")
        {
            fdl_pcl = pcl;
        }
    }

    Vec_fp_4d fp_by_fty_pcl_county_CVI(up_vec.size());
    Vec_fp_2d inf_COM(up_vec.size());
    if(infFarms.empty())
    {
        for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
        {
            Farm_type* fty = up_vec[fty_idx].getFarmType();
            std::vector<Prem_class*> pcl_vec = up_vec[fty_idx].getPremClasses();
            Vec_fp_3d fp_by_pcl_county(pcl_vec.size());
            for(int pcl_idx=0; pcl_idx<int(pcl_vec.size()); ++pcl_idx)
            {
                Vec_fp_2d fp_by_county(allCounties.size());
                for(size_t county_idx=0; county_idx<allCounties.size(); ++county_idx)
                {
                    County* county = allCounties[county_idx];
                    std::vector<Farm*> fvec;
                    for(Farm* f : county->get_premises(fty))
                    {
                        if(f->get_prem_class_idx() == pcl_idx)
                        {
                            if(f->getCOMIdentifier() == nullptr)
                            {

                                fvec.push_back(f);
                            }
                            else
                            {
                                inf_COM[fty_idx].push_back(f);
                            }
                        }
                    }
                    fp_by_county[county_idx] = fvec;

                }
                fp_by_pcl_county[pcl_idx].swap(fp_by_county);
            }
            fp_by_fty_pcl_county_CVI[fty_idx].swap(fp_by_pcl_county);
        }
    }
    else
    {
        for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
        {
            std::vector<Prem_class*> pcl_vec = up_vec[fty_idx].getPremClasses();
            fp_by_fty_pcl_county_CVI[fty_idx].resize(pcl_vec.size(), Vec_fp_2d(allCounties.size()));
        }

        for(Farm* f : infFarms)
        {
            if(f->getCOMIdentifier() == nullptr)
            {
                Prem_class* pcl = f->get_prem_class();
                int pcl_idx = pcl->idx;
                int county_idx = f->get_parent_county()->get_idx();
                if(pcl == mkt_pcl or pcl == fdl_pcl) //This is a market or feedlot and needs to send shipments to all farm types.
                {
                    for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
                    {
                        fp_by_fty_pcl_county_CVI[fty_idx][pcl_idx][county_idx].push_back(f);
                    }
                }
                else
                {
                    size_t fty_idx = f->get_farm_type()->get_index();
                    fp_by_fty_pcl_county_CVI[fty_idx][pcl_idx][county_idx].push_back(f);
                }
            }
            else
            {
                size_t fty_idx = f->get_farm_type()->get_index();
                inf_COM[fty_idx].push_back(f);
            }
        }
    }

    for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
    {
        USAMMv3_parameters& up = up_vec[fty_idx];
        Farm_type* fty = up.getFarmType();
        USAMMv3_period_parstruct& period_pars = up.get_parstruct_by_period_idx(time_period_idx);

        Vec_fp_3d& fp_by_pcl_county = fp_by_fty_pcl_county_CVI[fty_idx];
        std::vector<Prem_class*> pcl_vec = up.getPremClasses();
        for(size_t o_pcl_idx=0; o_pcl_idx<pcl_vec.size(); ++o_pcl_idx)
        {
            Prem_class* o_pcl = pcl_vec[o_pcl_idx];
            Vec_fp_2d& o_inf_vectors = fp_by_pcl_county[o_pcl_idx]; //List of counties and associated vectors of infected premises of origin prem class.
            for(size_t d_pcl_idx=0; d_pcl_idx<pcl_vec.size(); ++d_pcl_idx)
            {
                Prem_class* d_pcl = pcl_vec[d_pcl_idx];
                Vec_d_2d& cc_rates = period_pars.eval_cc_rates_CVI[o_pcl_idx][d_pcl_idx];
                Vec_d_1d& c_o_rates = period_pars.eval_c_o_rates_CVI[o_pcl_idx][d_pcl_idx];

                size_t n_o_size_bins = up.getSizeBins(fty, o_pcl).size(); //Used to avoid calling .size() a million times on the evaluated n-vectors further down.
                size_t n_d_size_bins = up.getSizeBins(fty, d_pcl).size();
                for(size_t o_county_idx=0; o_county_idx<o_inf_vectors.size(); ++o_county_idx)
                {
                    Vec_fp_1d& o_county_inf_vec = o_inf_vectors[o_county_idx];
                    if(o_county_inf_vec.empty())
                    {
                        continue;
                    }
                    double& c_o_rate = c_o_rates[o_county_idx];
                    if(c_o_rate < 0) //Not yet calculated.
                    {
                        c_o_rate = 0.0;
                        County* o_county = allCounties[o_county_idx];
                        State* o_state = o_county->get_parent_state();
                        int o_state_fips = o_state->get_code();
                        Vec_d_1d& dc_rates = cc_rates[o_county_idx];
                        const std::vector<double>& o_eval_N_vec = o_county->get_eval_N_vec_origin_CVI(fty_idx, o_pcl_idx, up, time_period_idx);
                        double ocov = o_county->get_county_ocov_weight(fty);
                        double outflow = period_pars.outflowCVI_vec[o_state_fips];
                        double outgoing_county_weight = 0.0;
                        for(size_t d_county_idx=0; d_county_idx<allCounties.size(); ++d_county_idx)
                        {
                            double& cc_rate = dc_rates[d_county_idx];
                            County* d_county = allCounties[d_county_idx];
                            State* d_state = d_county->get_parent_state();
                            const std::vector<double>& d_eval_N_vec = d_county->get_eval_N_vec_dest_CVI(fty_idx, d_pcl_idx, up, time_period_idx);
                            //Sum of outer product of the evaluated N-vectors (N prems in each bin size * evaluated weight of that bin size) is the county-county premsize rep.
                            double cc_premsize_rep = 0.0;
                            for(size_t i=0; i<n_o_size_bins; ++i)
                            {
                                double ow = o_eval_N_vec[i];
                                for(size_t j=0; j<n_d_size_bins; ++j)
                                {
                                    double dw = d_eval_N_vec[j];
                                    cc_premsize_rep += ow*dw;
                                }
                            }
                            //If within county and within the same prem class, adjust cc_premseize_rep for the rate of i sending to i.
                            if(o_county == d_county and o_pcl_idx == d_pcl_idx)
                            {
                                cc_premsize_rep -= d_county->get_within_county_adj_terms_CVI(fty_idx, d_pcl_idx, up, time_period_idx);
                            }

                            double& k_val_CVI = period_pars.county_eval_kernel_CVI[o_county_idx][d_county_idx];
                            if(k_val_CVI < 0.0)
                            {
                                k_val_CVI = up.kernel_fun(county_distance_matrix[o_county_idx][d_county_idx],
                                                          period_pars.aCVI_vec.at(o_state_fips),
                                                          period_pars.bCVI_vec.at(o_state_fips));
                            }

                            double c = period_pars.c_vec[pairing_function(o_pcl_idx, d_pcl_idx)];
                            double dcov = d_county->get_county_dcov_weight(fty);
                            double inflow = period_pars.inflowCVI_vec[d_state->get_code()];
                            cc_rate = cc_premsize_rep * c * ocov * dcov * outflow * inflow * k_val_CVI * delta_t;
                            c_o_rate += cc_rate;

                            // Also calculate the total destination prem weight per prem class for each county.
                            // This is used when calculating the premises' slaughter shipment rate and
                            // is needed when simulating bTB with slaughter surveillance.
                            double dcounty_d_prem_weight = 0.0;
                            for(size_t j=0; j<n_d_size_bins; ++j)
                            {
                                double dw = d_eval_N_vec[j];
                                dcounty_d_prem_weight += dw;
                            }
                            outgoing_county_weight += dcounty_d_prem_weight * c * ocov * dcov * outflow * inflow * k_val_CVI * delta_t;
                        }
                        o_county->set_outgoing_county_CVI_weight(fty_idx, o_pcl_idx, d_pcl_idx, outgoing_county_weight);
                    }

                    int n_c = draw_poisson(c_o_rate); //The total number of shipments out from the ocounty to ALL other counties.
                    if(n_c > 0)
                    {
                        County* o_county = allCounties[o_county_idx];

                        //Distribute the shipments among destination counties.
                        const Vec_d_1d& dcounty_weights = cc_rates[o_county_idx]; //These are the rates from ocounty to all other individual counties. We'll use them as weights to determine where the shipments go.
                        std::vector<unsigned int> n_cc_vec;
                        draw_multinomial(n_c, dcounty_weights, n_cc_vec);

                        //Make a vector of the weights of all the individual infected premises, and one additional element that represents all other premises (which we dont care about).
                        const std::vector<double>& ocounty_internal_prem_weights = o_county->get_USAMMv3_origin_prem_weights_CVI(fty_idx, o_pcl, time_period, up);
                        std::vector<double> ocounty_internal_inf_prem_weights(o_county_inf_vec.size()+1, 0.0);
                        ocounty_internal_inf_prem_weights[0] = o_county->get_USAMMv3_origin_prem_weight_sum_CVI(fty_idx, o_pcl, time_period, up); //This is the total sum of all origin prem weights...
                        for(size_t i=0; i<o_county_inf_vec.size(); ++i)
                        {
                            Farm* prem = o_county_inf_vec[i];
                            int prem_idx = prem->get_idx_in_county();
                            double w = ocounty_internal_prem_weights[prem_idx];
                            ocounty_internal_inf_prem_weights[0] -= w; //Subtract this individual infected prems' weight from the total representing all non-infected premises.
                            ocounty_internal_inf_prem_weights[i+1] = w;
                        }

                        for(size_t d_county_idx=0; d_county_idx<n_cc_vec.size(); ++d_county_idx)
                        {
                            unsigned int n_cc = n_cc_vec[d_county_idx];
                            if(n_cc > 0)
                            {
                                unsigned int n_processed_shipments = 0;
                                County* d_county = allCounties[d_county_idx];
                                const std::vector<double>& dcounty_internal_prem_weights = d_county->get_USAMMv3_dest_prem_weights_CVI(fty_idx, d_pcl, time_period, up);
                                const std::vector<Farm*>& dprem_ptrs = d_county->get_premises();
                                //Since the prem-prem weights only depend on prem size at this point (not even c) we can first sample origin premises,
                                //and then sample destination premises independently of each other for each shipment.
                                std::vector<unsigned int> o_prem_n_shipments;
                                draw_multinomial(n_cc, ocounty_internal_inf_prem_weights, o_prem_n_shipments);
                                for(size_t o_prem_idx=0; o_prem_idx<o_county_inf_vec.size(); ++o_prem_idx)
                                {
                                    size_t o_n = o_prem_n_shipments[o_prem_idx+1]; //Number of shipment originating from o_prem (+1 because element 0 represents all other prems).
                                    if(o_n > 0)
                                    {
                                        Farm* oprem = o_county_inf_vec[o_prem_idx];
                                        std::vector<unsigned int> d_prem_n_shipments;
                                        if(o_county == d_county and o_pcl == d_pcl)
                                        {
                                            //Remove the weight of the origin prem itself to prevent shipments to itself.
                                            std::vector<double> dcounty_internal_prem_weights_copy = dcounty_internal_prem_weights;
                                            dcounty_internal_prem_weights_copy[oprem->get_idx_in_county()] = 0.0;
                                            draw_multinomial(o_n, dcounty_internal_prem_weights_copy, d_prem_n_shipments);
                                        }
                                        else
                                        {
                                            draw_multinomial(o_n, dcounty_internal_prem_weights, d_prem_n_shipments);
                                        }

                                        for(size_t dprem_idx=0; dprem_idx<d_prem_n_shipments.size(); ++dprem_idx)
                                        {
                                            size_t d_n = d_prem_n_shipments[dprem_idx]; //Number of shipments going to dprem.
                                            if(d_n > 0)
                                            {
                                                Farm* dprem = dprem_ptrs[dprem_idx];
                                                if(oprem == dprem)
                                                {
                                                    std::cout << "Shipment origin == receiver." << std::endl;
                                                }
                                                //Origin and destination premises determined, now determine shipment size.
                                                for(size_t ship_idx=0; ship_idx<d_prem_n_shipments[dprem_idx]; ++ship_idx)
                                                {
                                                    size_t shipment_volume;
                                                    if(up.getUSAMMv3SpeciesCode() == "s")
                                                    {
                                                        shipment_volume = make_swine_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up);
                                                    }
                                                    else
                                                    {
                                                        shipment_volume = make_cattle_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up, mkt_pcl);
                                                    }

                                                    Shipment* s = new Shipment{int(timestep), // timestep of shipment
                                                                               day_of_year,
                                                                               oprem,
                                                                               dprem,
                                                                               o_county->get_id(),
                                                                               d_county->get_id(),
                                                                               fty->get_species(),
                                                                               static_cast<int>(shipment_volume),
                                                                               -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
                                                                               0, // ban (filled in filter_shipExposure if applicable)
                                                                               time_period,
                                                                               o_county->get_parent_state()->get_code(),
                                                                               o_county->get_parent_state()->get_id(),
                                                                               d_county->get_parent_state()->get_code(),
                                                                               d_county->get_parent_state()->get_id(),
                                                                               oprem->get_prem_class()->tag,
                                                                               oprem->get_size_allSpecies(),
                                                                               dprem->get_prem_class()->tag,
                                                                               dprem->get_size_allSpecies(),
                                                                               nullptr,
                                                                               "CVI"};
                                                    output.push_back(s);
                                                    ++n_processed_shipments;
                                                    ++n_ships_generated;
                                                }
                                            }
                                            if(n_processed_shipments >= n_cc)
                                            {
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // Now set the outgoing shipment rates of the individual premises that need
        // it later on when doing slaughter surveillance. (This cant be done from
        // any other location easily).
        for(Farm* f : infFarms)
        {
            int o_pcl_idx = f->get_prem_class_idx();
            size_t fty_idx = f->get_farm_type()->get_index();
            double prem_oweight = std::exp( f->get_USAMMv3_log_normed_binned_size()
                                  * period_pars.phiO_vec.at(o_pcl_idx) );
            double w = 0.0;
            for(size_t d_pcl_idx=0; d_pcl_idx<pcl_vec.size(); ++d_pcl_idx)
            {
                County* o_county = f->get_parent_county();
                double county_oweight =
                    o_county->get_outgoing_county_CVI_weight(fty_idx,
                                                             o_pcl_idx,
                                                             d_pcl_idx);
                w += prem_oweight * county_oweight;
            }
            f->set_latest_shipping_rate(w);
        }
    }



    Vec_d_2d evaluated_COM_prems_oweights(up_vec.size());
    Vec_d_2d evaluated_COM_prems_dweights(up_vec.size());
    for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
    {
        USAMMv3_parameters& up = up_vec[fty_idx];
        if(!up.get_COM_ids().empty())
        {
            USAMMv3_period_parstruct& period_pars = up.get_parstruct_by_period_idx(time_period_idx);
            std::vector<COM_prem_id>& COM_ids = up.get_COM_ids();

            Vec_d_2d& random_effects = period_pars.COM_random_effects.at(fty_idx);
            if(random_effects.empty())
            {
                random_effects.resize(COM_ids.size(), Vec_d_1d(COM_ids.size(), -std::numeric_limits<double>::max()));
            }
            Vec_d_1d& oweights_vec = evaluated_COM_prems_oweights[fty_idx];
            Vec_d_1d& dweights_vec = evaluated_COM_prems_dweights[fty_idx];

            oweights_vec.resize(COM_ids.size(), 0.0);
            dweights_vec.resize(COM_ids.size(), 0.0);
            for(size_t idx=0; idx<COM_ids.size(); ++idx)
            {
                COM_prem_id& COM_id = COM_ids[idx];
                Farm* farm = COM_id.get_farm();
                Prem_class* pcl = farm->get_prem_class();
                double normed_prem_size = farm->get_USAMMv3_binned_size() / up.get_avg_prem_size(pcl);
                double ophi = up.get_parstruct_by_period_idx(time_period_idx).phiO_vec.at(pcl->idx);
                double dphi = up.get_parstruct_by_period_idx(time_period_idx).phiD_vec.at(pcl->idx);
                oweights_vec[idx] = std::pow(normed_prem_size, ophi);
                dweights_vec[idx] = std::pow(normed_prem_size, dphi);
            }

            for(Farm* o_prem : inf_COM[fty_idx])
            {
                COM_prem_id* o_COM_id = o_prem->getCOMIdentifier();
                size_t o_COM_idx = o_COM_id->get_idx();
                Farm_type* o_fty = o_prem->get_farm_type();
                Prem_class* o_pcl = o_prem->get_prem_class();
                County* o_county = o_prem->get_parent_county();
                State* o_state = o_county->get_parent_state();

                int this_producer = o_COM_id->get_producer();
                Vec_d_1d& ostate_p_obs_COM_lookup_table = p_obs_COM_lookup_table[o_state->get_code()];

                //Calculate everything that only depends on the origin prem.
                double rho = period_pars.rho_vec.at(0);
                double oprem_size_w = evaluated_COM_prems_oweights[fty_idx][o_COM_idx];
                double outflow_COM = period_pars.outflowCOM_vec.at(o_state->get_code());
                double ocov = o_county->get_county_ocov_weight(o_fty);
                int o_county_idx = o_county->get_idx();
                double lambda_origin_part = rho * oprem_size_w * ocov * outflow_COM;
                Vec_d_1d& o_random_effects = random_effects[o_COM_idx];

                for(COM_prem_id* d_COM_id : up.get_COM_id_ptrs_by_producer(this_producer)) //Go through all other COM_ids that belong to the same producer.
                {
                    if(o_COM_id == d_COM_id)
                    {
                        continue; //No shipment to the prem itself.
                    }
                    size_t d_COM_idx = d_COM_id->get_idx();
                    Farm* dprem = d_COM_id->get_farm();
                    Prem_class* d_pcl = dprem->get_prem_class();
                    County* d_county = dprem->get_parent_county();
                    int d_county_idx = d_county->get_idx();
                    State* d_state = d_county->get_parent_state();
                    int d_state_fips = d_state->get_code();
                    double dprem_size_w = evaluated_COM_prems_dweights[fty_idx][d_COM_idx];
                    double dcov = d_county->get_county_dcov_weight(o_fty);
                    double inflow_COM = period_pars.inflowCOM_vec[d_state_fips];
                    double k_val_COM = period_pars.county_eval_kernel_COM[o_county_idx][d_county_idx];
                    if(k_val_COM < 0.0)
                    {
                        k_val_COM = up.kernel_fun(county_distance_matrix[o_county_idx][d_county_idx],
                                                  period_pars.aCOM_vec[o_state->get_code()],
                                                  period_pars.bCOM_vec[o_state->get_code()]);
                        period_pars.county_eval_kernel_COM[o_county_idx][d_county_idx] = k_val_COM;
                    }

                    //Random effect sampled here, only when it's actually needed in order to
                    //avoid sampling it for premises that aren't infected, and so wont't need it.
                    double p_obs_COM = ostate_p_obs_COM_lookup_table[d_state_fips];
                    double COM_rate_no_rnd_eff = lambda_origin_part * dprem_size_w *
                                                 dcov * inflow_COM * k_val_COM;

                    double& COM_rnd_effect_ij = o_random_effects[d_COM_idx];
                    if(COM_rnd_effect_ij == -std::numeric_limits<double>::max())
                    {
                        //Random effect has not been sampled for this pair yet.
                        int n_obs = o_COM_id->get_n_shipments_to(d_COM_id, time_period_idx);
                        double prod_alpha = period_pars.common_pars->alpha_map[this_producer];
                        double sample_shape = prod_alpha + n_obs;
                        double sample_rate = prod_alpha + COM_rate_no_rnd_eff * parameters->USAMM_temporal_n_timesteps[time_period_idx] * p_obs_COM; //Rho has been scaed down to be per day when reading posterior. Needs to be scaled up again for this purpose.
                        COM_rnd_effect_ij = gsl_ran_gamma(R, sample_shape, 1.0 / sample_rate);
                    }
                    double COM_rate = COM_rate_no_rnd_eff * COM_rnd_effect_ij;
                    int ij_n_shipments = gsl_ran_poisson(R, COM_rate);

                    for(int ship_idx=0; ship_idx<ij_n_shipments; ++ship_idx)
                    {
                        size_t shipment_volume = make_swine_shipment_volume(o_prem, dprem, o_pcl, d_pcl, time_period, up);
                        Shipment* s = new Shipment{int(timestep), // timestep of shipment
                                                   day_of_year,
                                                   o_prem,
                                                   dprem,
                                                   o_county->get_id(),
                                                   d_county->get_id(),
                                                   o_fty->get_species(),
                                                   static_cast<int>(shipment_volume),
                                                   -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
                                                   0, // ban (filled in filter_shipExposure if applicable)
                                                   time_period,
                                                   o_county->get_parent_state()->get_code(),
                                                   o_county->get_parent_state()->get_id(),
                                                   d_county->get_parent_state()->get_code(),
                                                   d_county->get_parent_state()->get_id(),
                                                   o_prem->get_prem_class()->tag,
                                                   o_prem->get_size_allSpecies(),
                                                   dprem->get_prem_class()->tag,
                                                   dprem->get_size_allSpecies(),
                                                   nullptr,
                                                   "COM"};
                        output.push_back(s);
                        ++n_ships_generated;
                        if(o_prem == dprem)
                        {
                            std::cout << "Shipment origin == receiver." << std::endl;
                        }
                    }
                } //End for d_COM_id.
            } //End for oprem.
        } //End if !COM_ids.empty()
    } //End for fty
}

//void Shipment_manager::makeShipmentsUSAMMv3(size_t timestep, size_t day_of_year, std::string time_period,
//                                            int time_period_idx, std::vector<Shipment*>& output,
//                                            std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec)
//{
//    Vec_d_2d evaluated_COM_prems_oweights(up_vec.size());
//    Vec_d_2d evaluated_COM_prems_dweights(up_vec.size());
//    for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
//    {
//        USAMMv3_parameters& up = up_vec[fty_idx];
//        Vec_d_1d& oweights_vec = evaluated_COM_prems_oweights[fty_idx];
//        Vec_d_1d& dweights_vec = evaluated_COM_prems_dweights[fty_idx];
//        std::vector<COM_prem_id>& COM_ids = up.get_COM_ids();
//        oweights_vec.resize(COM_ids.size(), 0.0);
//        dweights_vec.resize(COM_ids.size(), 0.0);
//        for(size_t idx=0; idx<COM_ids.size(); ++idx)
//        {
//            COM_prem_id& COM_id = COM_ids[idx];
//            Farm* farm = COM_id.get_farm();
//            Prem_class* pcl = farm->get_prem_class();
//            size_t fty_idx = farm->get_farm_type()->get_index();
//            USAMMv3_parameters& up = up_vec[fty_idx];
//            double normed_prem_size = farm->get_USAMMv3_binned_size() / up.get_avg_prem_size(pcl);
//            double ophi = up.get_parstruct_by_period_idx(time_period_idx).phiO_vec.at(pcl->idx);
//            double dphi = up.get_parstruct_by_period_idx(time_period_idx).phiD_vec.at(pcl->idx);
//            oweights_vec[idx] = std::pow(normed_prem_size, ophi);
//            dweights_vec[idx] = std::pow(normed_prem_size, dphi);
//        }
//    }
//
//    for(Farm* oprem : affected_farms)
//    {
//        int n_farm_shipments = 0;
//        double outgoing_shipment_rate = 0.0; //Keeps track of the sum of the rate of this prem sending to all other prems. For calculating the probability of sending to slaughter, used for the slaughter surveillance component.
//        Prem_class* o_pcl = oprem->get_prem_class();
////        int o_pcl_idx = o_pcl->idx;
//        std::vector<Farm_type*> o_ftypes;
//        if(o_pcl == mkt_pcl)
//        {
//            //If it's a market and cattle, generate shipments for both beef and dairy.
//            if(oprem->get_farm_type()->get_species() == "swine")
//            {
//               o_ftypes = { oprem->get_farm_type() };
//            }
//            else
//            {
//                for(auto& up : up_vec)
//                {
//                    if(up.getFarmType()->get_species() != "swine")
//                    {
//                        o_ftypes.push_back(up.getFarmType());
//                    }
//                }
//            }
//        }
//        else
//        {
//            o_ftypes = { oprem->get_farm_type() };
//        }
//
//        Farm_type* o_fty = oprem->get_farm_type();
//        County* o_county = oprem->get_parent_county();
//        State* o_state = o_county->get_parent_state();
//        USAMMv3_parameters& up = up_vec[o_fty->get_index()];
//        USAMMv3_period_parstruct& period_pars = up.get_parstruct_by_period_idx(time_period_idx);
//
//        if(oprem->getCOMIdentifier() == nullptr)
//        {
//            for(Farm_type* o_fty : o_ftypes)
//            {
//                int o_fty_idx = o_fty->get_index();
//                USAMMv3_parameters& up = up_vec[o_fty_idx];
//                std::vector<Prem_class*> ft_pclasses = up.getPremClasses();
//                County* o_county = oprem->get_parent_county();
//                Shipment_kernel kCVI(o_county->get_parent_state()->get_aCVI(o_fty),
//                                     o_county->get_parent_state()->get_bCVI(o_fty),
//                                     parameters->shipment_kernel, true);
//
//                double normed_opremsize = oprem->get_USAMMv3_binned_size() / up.get_avg_prem_size(o_pcl);
//                double ophi = up.get_phi_O(o_pcl, time_period);
//                double oprem_size_w = std::pow(normed_opremsize, ophi);
//                double outflow_CVI = o_county->get_parent_state()->get_outflowCVI(o_fty);
//                double ocov = o_county->get_county_ocov_weight(o_fty);
//
//                for(size_t d_county_idx=0; d_county_idx<allCounties.size(); ++d_county_idx)
//                {
//                    County* d_county = allCounties[d_county_idx];
//                    for(Prem_class* d_pcl : d_county->get_prem_classes_by_type_idx(o_fty_idx)) //Only go through the pclasses that are present in the receiving county.
//                    {
//                        double d_county_premsize_rep = d_county->get_dest_premsize_rep_CVI(o_fty, o_pcl, d_pcl, up, time_period, time_period_idx); //Includes all prems in dcounty, their size, dphi and c.
//
//
//                        if(d_county_premsize_rep == 0)
//                        {
//                            continue; //I'm not entirely sure this will ever happen, but this check is not done very often so I'll leave it.
//                        }
//
//                        if(o_county == d_county and o_pcl == d_pcl)
//                        {
//                            //Within county shipment needs to correct for the fact that a prem can't ship to itself.
//                            d_county_premsize_rep -= std::pow(normed_opremsize, up.get_phi_D(o_pcl, time_period)) *
//                                                     up.get_c(o_pcl, d_pcl, time_period);
//                        }
//
//                        double inflow_CVI = d_county->get_parent_state()->get_inflowCVI(o_fty);
//                        double dcov = d_county->get_county_dcov_weight(o_fty);
//
//                        double k_val_CVI = kernel_vals_CVI[o_county->get_idx()][d_county->get_idx()];
//                        if(k_val_CVI < 0.0)
//                        {
//                            k_val_CVI = kCVI.kernel(o_county, d_county);
//                            kernel_vals_CVI[o_county->get_idx()][d_county->get_idx()] = k_val_CVI;
//                        }
//
//                        double d_county_receiving_rate = d_county_premsize_rep *
//                                                         oprem_size_w *
//                                                         outflow_CVI *
//                                                         inflow_CVI *
//                                                         ocov *
//                                                         dcov *
//                                                         k_val_CVI;
//                        int n_county_shipments = draw_poisson(d_county_receiving_rate);
//                        if(d_county_receiving_rate < 0.0)
//                        {
//                            std::cout << "Negative shipping rate." << std::endl;
//                        }
//                        outgoing_shipment_rate += d_county_receiving_rate;
//                        if(n_county_shipments > 0)
//                        {
//                            int n_processed_shipments = 0;
//                            std::vector<double> dcounty_internal_prem_weights = d_county->get_USAMMv3_dest_prem_weights_CVI(o_fty, d_pcl, time_period, up_vec);
//                            std::vector<unsigned int> dprem_n_ships;
//                            if(o_county == d_county and o_pcl == d_pcl)
//                            {
//                                //Remove the weight of the origin prem itself to prevent shipments to itself.
//                                std::vector<double> dcounty_internal_prem_weights_copy = dcounty_internal_prem_weights;
//                                dcounty_internal_prem_weights_copy[oprem->get_idx_in_county()] = 0.0;
//                                draw_multinomial(n_county_shipments, dcounty_internal_prem_weights_copy, dprem_n_ships);
//                            }
//                            else
//                            {
//                                draw_multinomial(n_county_shipments, dcounty_internal_prem_weights, dprem_n_ships);
//                            }
//
//                            std::vector<Farm*> dprem_ptrs = d_county->get_premises();
//                            for(size_t dprem_idx=0; dprem_idx<dprem_n_ships.size(); ++dprem_idx)
//                            {
//                                if(dprem_n_ships[dprem_idx] > 0)
//                                {
//                                    Farm* dprem = dprem_ptrs[dprem_idx];
//                                    //Origin and destination premises determined, now determine shipment size.
//                                    for(size_t ship_idx=0; ship_idx<dprem_n_ships[dprem_idx]; ++ship_idx)
//                                    {
//                                        size_t shipment_volume;
//                                        if(up.getUSAMMv3SpeciesCode() == "s")
//                                        {
//                                            shipment_volume = make_swine_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up);
//                                        }
//                                        else
//                                        {
//                                            shipment_volume = make_cattle_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up, mkt_pcl);
//                                        }
//
//                                        Shipment* s = new Shipment{int(timestep), // timestep of shipment
//                                                                   day_of_year,
//                                                                   oprem,
//                                                                   dprem,
//                                                                   o_county->get_id(),
//                                                                   d_county->get_id(),
//                                                                   o_fty->get_species(),
//                                                                   static_cast<int>(shipment_volume),
//                                                                   -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
//                                                                   0, // ban (filled in filter_shipExposure if applicable)
//                                                                   time_period,
//                                                                   o_county->get_parent_state()->get_code(),
//                                                                   o_county->get_parent_state()->get_id(),
//                                                                   d_county->get_parent_state()->get_code(),
//                                                                   d_county->get_parent_state()->get_id(),
//                                                                   oprem->get_prem_class()->tag,
//                                                                   oprem->get_size_allSpecies(),
//                                                                   dprem->get_prem_class()->tag,
//                                                                   dprem->get_size_allSpecies(),
//                                                                   nullptr,
//                                                                   "CVI"};
//                                        output.push_back(s);
//                                        ++n_processed_shipments;
//                                        ++n_ships_generated;
//                                        ++n_farm_shipments;
//                                        if(oprem == dprem)
//                                        {
//                                            std::cout << "Shipment origin == receiver." << std::endl;
//                                        }
//                                    }
//                                }
//                                if(n_processed_shipments >= n_county_shipments)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//            oprem->set_latest_shipping_rate(outgoing_shipment_rate);
//        }
//        else
//        {
//            //It's a COM premises.
//            COM_prem_id* this_COM_id = oprem->getCOMIdentifier();
//            int this_producer = this_COM_id->get_producer();
//            Vec_d_1d& ostate_p_obs_COM_lookup_table = p_obs_COM_lookup_table[o_state->get_code()];
//
//            //Calculate everything that only depends on the origin prem.
//            double rho = period_pars.rho_vec.at(0);
//            double oprem_size_w = evaluated_COM_prems_oweights[o_fty->get_index()][this_COM_id->get_idx()];
//            double outflow_COM = period_pars.outflowCOM_vec.at(o_state->get_code());
//            double ocov = o_county->get_county_ocov_weight(o_fty);
//            double lambda_origin_part = rho * oprem_size_w * ocov * outflow_COM;
//
//            for(COM_prem_id* other_COM_id : up.get_COM_id_ptrs_by_producer(this_producer)) //Go through all other COM_ids that belong to the same producer.
//            {
//                if(this_COM_id == other_COM_id)
//                {
//                    continue; //No shipment to the prem itself.
//                }
//                Farm* dprem = other_COM_id->get_farm();
//                Prem_class* d_pcl = dprem->get_prem_class();
//                County* d_county = dprem->get_parent_county();
//                State* d_state = d_county->get_parent_state();
//                double dprem_size_w = evaluated_COM_prems_dweights[o_fty->get_index()][other_COM_id->get_idx()];
//                double dcov = d_county->get_county_dcov_weight(o_fty);
//                double inflow_COM = period_pars.inflowCOM_vec.at(d_state->get_code());
//                double k_val_COM = period_pars.county_eval_kernel_COM[o_county->get_idx()][d_county->get_idx()];
//                if(k_val_COM < 0.0)
//                {
//                    k_val_COM = up.kernel_fun(county_distance_matrix[o_county->get_idx()][d_county->get_idx()],
//                                              period_pars.aCOM_vec.at(o_state->get_code()),
//                                              period_pars.bCOM_vec.at(o_state->get_code()));
//                    period_pars.county_eval_kernel_COM[o_county->get_idx()][d_county->get_idx()] = k_val_COM;
//                }
//
//                //Random effect sampled here, only when it's actually needed in order to
//                //avoid sampling it for premises that aren't infected, and so wont't need it.
//                double p_obs_COM = ostate_p_obs_COM_lookup_table.at(d_state->get_code());
//                double COM_rate_no_rnd_eff = lambda_origin_part * dprem_size_w *
//                                             dcov * inflow_COM * k_val_COM;
//
//                double COM_rnd_effect_ij = this_COM_id->get_random_effect(other_COM_id->get_idx(), time_period_idx);
//                if(COM_rnd_effect_ij == -std::numeric_limits<double>::max())
//                {
//                    //Random effect has not been sampled for this pair yet.
//                    int n_obs = this_COM_id->get_n_shipments_to(other_COM_id, time_period_idx);
//                    double prod_alpha = period_pars.common_pars->alpha_map.at(this_producer);
//                    double sample_shape = prod_alpha + n_obs;
//                    double sample_rate = prod_alpha + COM_rate_no_rnd_eff * parameters->USAMM_temporal_n_timesteps[time_period_idx] * p_obs_COM; //Rho has been scaed down to be per day when reading posterior. Needs to be scaled up again for this purpose.
//                    COM_rnd_effect_ij = gsl_ran_gamma(R, sample_shape, 1.0 / sample_rate);
//                    this_COM_id->set_random_effect(other_COM_id->get_idx(), time_period_idx, COM_rnd_effect_ij);
//                }
//                double COM_rate = COM_rate_no_rnd_eff * COM_rnd_effect_ij;
//                int ij_n_shipments = gsl_ran_poisson(R, COM_rate);
//
//                for(int ship_idx=0; ship_idx<ij_n_shipments; ++ship_idx)
//                {
//                    size_t shipment_volume = make_swine_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up);
//                    Shipment* s = new Shipment{int(timestep), // timestep of shipment
//                                               day_of_year,
//                                               oprem,
//                                               dprem,
//                                               o_county->get_id(),
//                                               d_county->get_id(),
//                                               o_fty->get_species(),
//                                               static_cast<int>(shipment_volume),
//                                               -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
//                                               0, // ban (filled in filter_shipExposure if applicable)
//                                               time_period,
//                                               o_county->get_parent_state()->get_code(),
//                                               o_county->get_parent_state()->get_id(),
//                                               d_county->get_parent_state()->get_code(),
//                                               d_county->get_parent_state()->get_id(),
//                                               oprem->get_prem_class()->tag,
//                                               oprem->get_size_allSpecies(),
//                                               dprem->get_prem_class()->tag,
//                                               dprem->get_size_allSpecies(),
//                                               nullptr,
//                                               "COM"};
//                    output.push_back(s);
//                    ++n_ships_generated;
//                    if(oprem == dprem)
//                    {
//                        std::cout << "Shipment origin == receiver." << std::endl;
//                    }
//                }
//            }
//        }
//    } //End for farm.

//    std::clock_t t_start = std::clock();
//
//    std::clock_t t_end = std::clock();
//    double t_tot = 1000.0 * (t_end - t_start) / CLOCKS_PER_SEC;
//    std::cout << "t_tot: " << t_tot << "ms." << std::endl;
//    exit(0);
//}


//void Shipment_manager::makeShipmentsUSAMMv3_new(size_t timestep, size_t day_of_year, std::string time_period,
//                                            int time_period_idx, std::vector<Shipment*>& output,
//                                            std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec)
//{
//    int n_ships_generated = 0;
//    std::vector<std::vector<double>> kernel_vals_CVI;
//    kernel_vals_CVI.resize(allCounties.size(), std::vector<double>(allCounties.size(), -1.0));
//
//    std::vector<std::vector<std::vector<Prem_class*>>> prem_classes_by_fty_and_county;
//    prem_classes_by_fty_and_county.resize(up_vec.size(), std::vector<std::vector<Prem_class*>>(allCounties.size()));
//    for(size_t fty_idx=0; fty_idx<up_vec.size(); ++fty_idx)
//    {
//        for(size_t county_idx=0; county_idx<allCounties.size(); ++county_idx)
//        {
//            County* c = allCounties[county_idx];
//            prem_classes_by_fty_and_county[fty_idx][county_idx] = c->get_prem_classes_by_type_idx(fty_idx);
//        }
//    }
//
//    Vec_d_4d c_premsize_reps(up_vec.size());
//    for(USAMMv3_parameters& up : up_vec)
//    {
//        Farm_type* fty = up.getFarmType();
//        size_t n_pcl = up.getPremClasses().size();
//        Vec_d_3d v(n_pcl, Vec_d_2d(n_pcl, Vec_d_1d(allCounties.size())));
//        for(size_t o_pcl_idx=0; o_pcl_idx<n_pcl; ++o_pcl_idx)
//        {
//            Prem_class* o_pcl = up.getPremClasses()[o_pcl_idx];
//            for(size_t d_pcl_idx=0; d_pcl_idx<n_pcl; ++d_pcl_idx)
//            {
//                Prem_class* d_pcl = up.getPremClasses()[d_pcl_idx];
//                for(size_t c_idx=0; c_idx<allCounties.size(); ++c_idx)
//                {
//                    County* c = allCounties[c_idx];
//                    v[o_pcl_idx][d_pcl_idx][c_idx] = c->get_dest_premsize_rep_CVI(fty, o_pcl, d_pcl, up,
//                                                                                 time_period, time_period_idx);
//
//                }
//            }
//        }
//        c_premsize_reps[fty->get_index()] = v;
//    }
//
//    //Select what farms will be involved in the generation of shipments.
//    //Running with an empty infFarms is a signal to generating a full network of shipments,
//    //so then we use all farms.
//    //This function does have anything to do with btb, as btb is always simulated with USAMMv3
//    std::vector<Farm*> affected_farms;
//    affected_farms.reserve(1000000);
//    if(infFarms.empty())
//    {
//        for(County* c : allCounties)
//        {
//            std::vector<Farm*> c_farms = c->get_premises();
//            affected_farms.insert(affected_farms.end(), c_farms.begin(), c_farms.end());
//        }
//    }
//    else
//    {
//        affected_farms.assign(infFarms.begin(), infFarms.end());
//    }
//
//    //Determine the index of the market prem class (just for convenience/speed).
//    Prem_class* mkt_pcl = nullptr;
//    std::vector<Prem_class*> pcl_vec = up_vec.begin()->getPremClasses();
//    for(Prem_class* pcl : pcl_vec)
//    {
//        if(pcl->tag == "Mkt")
//        {
//            mkt_pcl = pcl;
//            break;
//        }
//    }
//
//
//    //This function is only for making shipments during a disease simulation.
//    for(Farm* oprem : affected_farms)
//    {
//        int n_farm_shipments = 0;
//        double outgoing_shipment_rate = 0.0; //Keeps track of the sum of the rate of this prem sending to all other prems. For calculating the probability of sending to slaughter, used for the slaughter surveillance component.
//        Prem_class* o_pcl = oprem->get_prem_class();
//        std::vector<Farm_type*> o_ftypes;
//        if(o_pcl == mkt_pcl)
//        {
//            //If it's a market and cattle, generate shipments for both beef and dairy.
//            if(oprem->get_farm_type()->get_species() == "swine")
//            {
//               o_ftypes = { oprem->get_farm_type() };
//            }
//            else
//            {
//                for(auto& up : up_vec)
//                {
//                    if(up.getFarmType()->get_species() != "swine")
//                    {
//                        o_ftypes.push_back(up.getFarmType());
//                    }
//                }
//            }
//
//        }
//        else
//        {
//            o_ftypes = { oprem->get_farm_type() };
//        }
//
//        if(oprem->getCOMIdentifier() == nullptr)
//        {
//            for(Farm_type* o_fty : o_ftypes)
//            {
//                int o_fty_idx = o_fty->get_index();
//                USAMMv3_parameters& up = up_vec[o_fty_idx];
//                std::vector<Prem_class*> ft_pclasses = up.getPremClasses();
//                County* o_county = oprem->get_parent_county();
//                Shipment_kernel kCVI(o_county->get_parent_state()->get_aCVI(o_fty),
//                                     o_county->get_parent_state()->get_bCVI(o_fty),
//                                     parameters->shipment_kernel, true);
//
//                double normed_opremsize = oprem->get_USAMMv3_binned_size() / up.get_avg_prem_size(o_pcl);
//                double ophi = up.get_phi_O(o_pcl, time_period);
//                double oprem_size_w = std::pow(normed_opremsize, ophi);
//                double outflow_CVI = o_county->get_parent_state()->get_outflowCVI(o_fty);
//                double ocov = o_county->get_county_ocov_weight(o_fty);
//
//                std::vector<std::vector<Prem_class*>>& fty_c_pcl_vec = prem_classes_by_fty_and_county[o_fty_idx];
//                for(size_t d_county_idx=0; d_county_idx<allCounties.size(); ++d_county_idx)
//                {
//                    County* d_county = allCounties[d_county_idx];
////                    for(Prem_class* d_pcl : d_county->get_prem_classes_by_type_idx(o_fty_idx)) //Only go through the pclasses that are present in the receiving county.
//                    auto& county_pcl_vec = fty_c_pcl_vec[d_county_idx];
//                    for(size_t c_pcl_idx=0; c_pcl_idx<county_pcl_vec.size(); ++c_pcl_idx) //Only go through the pclasses that are present in the receiving county.
//                    {
//                        Prem_class* d_pcl = county_pcl_vec[c_pcl_idx];
////                        double d_county_premsize_rep = d_county->get_dest_premsize_rep_CVI(o_fty, o_pcl, d_pcl, up, time_period, time_period_idx); //Includes all prems in dcounty, their size, dphi and c.
//                        double d_county_premsize_rep = c_premsize_reps[o_fty_idx][o_pcl->idx][d_pcl->idx][d_county_idx];
//
//                        if(d_county_premsize_rep == 0)
//                        {
//                            continue; //I'm not entirely sure this will ever happen, but this check is not done very often so I'll leave it.
//                        }
//
//                        if(o_county == d_county and o_pcl == d_pcl)
//                        {
//                            //Within county shipment needs to correct for the fact that a prem can't ship to itself.
//                            d_county_premsize_rep -= std::pow(normed_opremsize, up.get_phi_D(o_pcl, time_period)) *
//                                                     up.get_c(o_pcl, d_pcl, time_period);
//                        }
//
//                        double inflow_CVI = d_county->get_parent_state()->get_inflowCVI(o_fty);
//                        double dcov = d_county->get_county_dcov_weight(o_fty);
//
//                        double k_val_CVI = kernel_vals_CVI[o_county->get_idx()][d_county->get_idx()];
//                        if(k_val_CVI < 0.0)
//                        {
//                            k_val_CVI = kCVI.kernel(o_county, d_county);
//                            kernel_vals_CVI[o_county->get_idx()][d_county->get_idx()] = k_val_CVI;
//                        }
//
//                        double d_county_receiving_rate = d_county_premsize_rep *
//                                                         oprem_size_w *
//                                                         outflow_CVI *
//                                                         inflow_CVI *
//                                                         ocov *
//                                                         dcov *
//                                                         k_val_CVI;
//                        int n_county_shipments = draw_poisson(d_county_receiving_rate);
//                        if(d_county_receiving_rate < 0.0)
//                        {
//                            std::cout << "Negative shipping rate." << std::endl;
//                        }
//                        outgoing_shipment_rate += d_county_receiving_rate;
//                        if(n_county_shipments > 0)
//                        {
//                            int n_processed_shipments = 0;
//                            std::vector<double> dcounty_internal_prem_weights = d_county->get_USAMMv3_dest_prem_weights_CVI(o_fty, d_pcl, time_period, up_vec);
//                            std::vector<unsigned int> dprem_n_ships;
//                            if(o_county == d_county and o_pcl == d_pcl)
//                            {
//                                //Remove the weight of the origin prem itself to prevent shipments to itself.
//                                std::vector<double> dcounty_internal_prem_weights_copy = dcounty_internal_prem_weights;
//                                dcounty_internal_prem_weights_copy[oprem->get_idx_in_county()] = 0.0;
//                                draw_multinomial(n_county_shipments, dcounty_internal_prem_weights_copy, dprem_n_ships);
//                            }
//                            else
//                            {
//                                draw_multinomial(n_county_shipments, dcounty_internal_prem_weights, dprem_n_ships);
//                            }
//
//                            std::vector<Farm*> dprem_ptrs = d_county->get_premises();
//                            for(size_t dprem_idx=0; dprem_idx<dprem_n_ships.size(); ++dprem_idx)
//                            {
//                                if(dprem_n_ships[dprem_idx] > 0)
//                                {
//                                    Farm* dprem = dprem_ptrs[dprem_idx];
//                                    //Origin and destination premises determined, now determine shipment size.
//                                    for(size_t ship_idx=0; ship_idx<dprem_n_ships[dprem_idx]; ++ship_idx)
//                                    {
//                                        size_t shipment_volume;
//                                        if(up.getUSAMMv3SpeciesCode() == "s")
//                                        {
//                                            shipment_volume = make_swine_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up);
//                                        }
//                                        else
//                                        {
//                                            shipment_volume = make_cattle_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up, mkt_pcl);
//                                        }
//
//                                        Shipment* s = new Shipment{int(timestep), // timestep of shipment
//                                                                   day_of_year,
//                                                                   oprem,
//                                                                   dprem,
//                                                                   o_county->get_id(),
//                                                                   d_county->get_id(),
//                                                                   o_fty->get_species(),
//                                                                   static_cast<int>(shipment_volume),
//                                                                   -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
//                                                                   0, // ban (filled in filter_shipExposure if applicable)
//                                                                   time_period,
//                                                                   o_county->get_parent_state()->get_code(),
//                                                                   o_county->get_parent_state()->get_id(),
//                                                                   d_county->get_parent_state()->get_code(),
//                                                                   d_county->get_parent_state()->get_id(),
//                                                                   oprem->get_prem_class()->tag,
//                                                                   oprem->get_size_allSpecies(),
//                                                                   dprem->get_prem_class()->tag,
//                                                                   dprem->get_size_allSpecies(),
//                                                                   nullptr,
//                                                                   "CVI"};
//                                        output.push_back(s);
//                                        ++n_processed_shipments;
//                                        ++n_ships_generated;
//                                        ++n_farm_shipments;
//                                        if(oprem == dprem)
//                                        {
//                                            std::cout << "Shipment origin == receiver." << std::endl;
//                                        }
//                                    }
//                                }
//                                if(n_processed_shipments >= n_county_shipments)
//                                {
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//            oprem->set_latest_shipping_rate(outgoing_shipment_rate);
//        }
//        else
//        {
//            //It's a COM premises.
//            Farm_type* o_fty = oprem->get_farm_type();
//            USAMMv3_parameters& up = up_vec[o_fty->get_index()];
//            COM_prem_id* this_COM_id = oprem->getCOMIdentifier();
//            int this_producer = this_COM_id->get_producer();
//            for(COM_prem_id& other_COM_id : up.get_COM_ids())
//            {
//                if(this_COM_id == &other_COM_id)
//                {
//                    continue; //No shipment to the prem itself.
//                }
//                Farm* dprem = other_COM_id.get_farm();
//                County* o_county = oprem->get_parent_county();
//                County* d_county = dprem->get_parent_county();
//                if(this_producer == other_COM_id.get_producer() and
//                   o_county != d_county) //Same producer and different counties, so shipments can be sent between these.
//                {
//                    Prem_class* d_pcl = dprem->get_prem_class();
//                    State* o_state = o_county->get_parent_state();
//                    State* d_state = d_county->get_parent_state();
//
//                    double c = up.get_c(o_pcl, d_pcl, time_period);
//                    double rho = up.get_rho(o_pcl, d_pcl, time_period);
//                    double normed_opremsize = oprem->get_USAMMv3_binned_size() / up.get_avg_prem_size(o_pcl);
//                    double ophi = up.get_phi_O(o_pcl, time_period);
//                    double oprem_size_w = std::pow(normed_opremsize, ophi);
//                    double normed_dpremsize = dprem->get_USAMMv3_binned_size() / up.get_avg_prem_size(d_pcl);
//                    double dphi = up.get_phi_D(d_pcl, time_period);
//                    double dprem_size_w = std::pow(normed_dpremsize, dphi);
//                    double outflow_COM = o_county->get_parent_state()->get_outflowCOM(o_fty);
//                    double ocov = o_county->get_county_ocov_weight(o_fty);
//                    double dcov = d_county->get_county_dcov_weight(o_fty);
//                    double inflow_COM = d_county->get_parent_state()->get_inflowCOM(o_fty);
//                    Shipment_kernel kCOM(o_county->get_parent_state()->get_aCOM(o_fty),
//                                         o_county->get_parent_state()->get_bCOM(o_fty),
//                                         parameters->shipment_kernel, true);
//                    double k_val_COM = kCOM.kernel(o_county, d_county);
//
//                    //Random effect sampled here, only when it's actually needed in order to
//                    //avoid sampling it for premises that aren't infected, and so wont't need it.
//
//                    double p_obs_COM = p_obs_COM_lookup_table[o_state->get_code()][d_state->get_code()];
//                    double COM_rate_no_rnd_eff = c * rho * oprem_size_w * dprem_size_w *
//                                                 ocov * dcov * outflow_COM * inflow_COM * k_val_COM;
//
//                    double COM_rnd_effect_ij = this_COM_id->get_random_effect(other_COM_id.get_idx(), time_period_idx);
//                    if(COM_rnd_effect_ij == -std::numeric_limits<double>::max())
//                    {
//                        //Random effect has not been sampled for this pair yet.
//                        int n_obs = this_COM_id->get_n_shipments_to(&other_COM_id, time_period_idx);
//                        double prod_alpha = up.get_producer_alpha(this_producer);
//                        double sample_alpha = prod_alpha + n_obs;
//                        double sample_beta = prod_alpha + COM_rate_no_rnd_eff * p_obs_COM;
//                        COM_rnd_effect_ij = gsl_ran_gamma(R, sample_alpha, 1.0 / sample_beta);
//                    }
//                    double COM_rate = COM_rate_no_rnd_eff * COM_rnd_effect_ij;
//                    int ij_n_shipments = gsl_ran_poisson(R, COM_rate);
//
//                    for(int ship_idx=0; ship_idx<ij_n_shipments; ++ship_idx)
//                    {
//                        size_t shipment_volume = make_swine_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up);
//                        Shipment* s = new Shipment{int(timestep), // timestep of shipment
//                                                   day_of_year,
//                                                   oprem,
//                                                   dprem,
//                                                   o_county->get_id(),
//                                                   d_county->get_id(),
//                                                   o_fty->get_species(),
//                                                   static_cast<int>(shipment_volume),
//                                                   -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
//                                                   0, // ban (filled in filter_shipExposure if applicable)
//                                                   time_period,
//                                                   o_county->get_parent_state()->get_code(),
//                                                   o_county->get_parent_state()->get_id(),
//                                                   d_county->get_parent_state()->get_code(),
//                                                   d_county->get_parent_state()->get_id(),
//                                                   oprem->get_prem_class()->tag,
//                                                   oprem->get_size_allSpecies(),
//                                                   dprem->get_prem_class()->tag,
//                                                   dprem->get_size_allSpecies(),
//                                                   nullptr,
//                                                   "COM"};
//                        output.push_back(s);
//                        ++n_ships_generated;
//                        if(oprem == dprem)
//                        {
//                            std::cout << "Shipment origin == receiver." << std::endl;
//                        }
//                    }
//                }
//            }
//        }
//    } //End for farm.
//}

//void Shipment_manager::makeNetworkUSAMMv3(std::vector<std::string>out_fnames, Grid_manager& G)
//{
//    std::vector<std::ofstream*> f_vec;
//    for(std::string fname : out_fnames)
//    {
//        std::ofstream* f = new std::ofstream(fname + ".network");
//        if(f->is_open())
//        {
//            //Write the header to output file.
//            *f << "index\toCountyId\tdCountyId\tdayOfYear\tvolume\tunused\tperiod\toStateAbbr\toStateId\tdStateAbbr\tdStateId\toPremType\toPremSize\tdPremType\tdPremSize" <<
//                  //"\tdistance\tunnormProb\tdState_inflow\tdState_n\tdCounty_n\tdCounty_weight" <<
//                  std::endl;
//            f_vec.push_back(f);
//        }
//        else
//        {
//            std::cout << "Failed to open network generation output file: " << fname
//                      << ". Exiting..." << std::endl;
//            exit(EXIT_FAILURE);
//        }
//        delete f;
//    }
//
//
//    std::cout << "Generating shipments..." << std::endl;
//    std::vector<Shipment*> shipments;
//    shipments.reserve(1000000);
//    size_t shipment_counter = 1;
//    for(int day_i = 1; day_i < parameters->timesteps+1; day_i++)
//    {
//        int day_of_year = get_day_of_year(day_i, parameters->start_day);
//        G.updateShippingParameters(day_i, day_of_year);
//        std::string time_period = G.get_time_period();
//        size_t rem_days = G.get_rem_days_of_period();
//        std::vector<Farm*> dummy_vec = {};
//        std::cout << "Day " << day_i << std::endl;
//        makeShipmentsUSAMMv3(day_i, day_of_year, time_period, shipments, dummy_vec,
//                             G.get_usammv3_parameters_map());
//        std::cout << "Day " << day_i << " complete." << std::endl;
//        for(Shipment* s : shipments)
//        {
//            Farm_type* ft = G.get_farm_type_by_name(s->species);
//            size_t ft_index = ft->get_index();
//            *(f_vec.at(ft_index)) << shipment_counter << "\t"
//                        << s->origFIPS << "\t"
//                        << s->destFIPS << "\t"
//                        << s->day_of_year << "\t"
//                        << s->volume << "\t"
//                        << 1 << "\t"
//                        << time_period << "\t"
//                        << s->origState_abbrev << "\t"
//                        << s->origState_id << "\t"
//                        << s->destState_abbrev << "\t"
//                        << s->destState_id << "\t"
//                        << s->originIndType << "\t"
//                        << s->originSize << "\t"
//                        << s->destIndType << "\t"
//                        << s->destSize
////                        << 0.0 << "\t" //Distance
////                        << 0.0 << "\t" //Probability
////                        << 0.0 << "\t" //Dest state inflow
////                        << 0 << "\t" //Dest state n farms
////                        << 0 << "\t" //Dest county n farms
////                        << 0.0 << "\t" //Dest county inflow (weight).
////                        << s->origID << "\t" //Origin farm.
////                        << s->destID << "\t" //Destination farm.
//                        << std::endl;
//            shipment_counter += 1;
//        }
//        shipments.clear();
//    }
//    for(size_t i = 0; i > f_vec.size(); i++)
//    {
//        f_vec[i]->close();
//        delete f_vec[i];
//    }
//
//    std::vector<Farm_type*> ft_vec = G.get_farm_types();
//    for(size_t i = 0; i < ft_vec.size(); i++)
//    {
//        std::ofstream f(out_fnames.at(i) + ".gen");
//        if(f.is_open())
//        {
//            f << G.get_generation_string(ft_vec.at(i));
//            f.close();
//        }
//    }
//    std::cout << "...done." << std::endl;


    /////////////////////////////////////////////////////////////


//    std::cout << timestep << ", " << day_of_year << std::endl;
//    std::vector<unsigned int> rec_county_outcomes; //To store the number of shipments to each county. Declared here so it can be reused.
//    output.clear();
//    rec_county_outcomes.reserve(allCounties.size());
//
//    //Determine the index of the market prem class (just for convenience/speed).
//    Prem_class* mkt_pcl = nullptr;
//    std::vector<Prem_class*> pcl_vec = up_map.begin()->second.getPremClasses();
//    for(Prem_class* pcl : pcl_vec)
//    {
//        if(pcl->tag == "Mkt")
//        {
//            mkt_pcl = pcl;
//            break;
//        }
//    }
//
//    for(auto& ft_up_pair : up_map)
//    {
//        Farm_type* ft = ft_up_pair.first;
//        USAMMv3_parameters& up = up_map[ft];
//        std::vector<Prem_class*> ft_prclasses = up.getPremClasses();
//        for(Prem_class* o_pcl : ft_prclasses)
//        {
//            for(Prem_class* d_pcl : ft_prclasses)
//            {
//                for(size_t o_county_idx=0; o_county_idx<allCounties.size(); ++o_county_idx)
//                {
//                    County* o_county = allCounties[o_county_idx];
//                    int fips = o_county->get_fips_code();
//                    //This is the total rate of shipments from o_pcl to d_pcl leaving this county for any other county (or itself).
//                    double daily_ocounty_shipment_rate = o_county->get_o_CVI_shipment_rate(ft, o_pcl, d_pcl, time_period, up_map);
//                    size_t n_shipments = gsl_ran_poisson(R, daily_ocounty_shipment_rate);
//
//                    if(n_shipments > 0)
//                    {
//                        //Determine destination counties.
//                        std::vector<double>& internal_prem_weights_o = o_county->get_USAMMv3_origin_prem_weights_CVI(ft, o_pcl, time_period, up_map);
//                        std::vector<Farm*> ocounty_premises = o_county->get_premises();
//                        std::vector<unsigned int> oprem_outcome; //To store the number of shipments from different premises.
//                        oprem_outcome.reserve(internal_prem_weights_o.size());
//                        //These are the shipment rates from this county to each and every other county (incl. itself)
//                        //for shipments of the specific combination of prem classes. They are rates, but the probabilities
//                        //of the respective counties being receivers of the shipments are proportional to the rates. So
//                        //see them as unnormalized probabilities (weights). The multinomial function normalizes them internally.
//                        const std::vector<double>& d_lambdas = o_county->get_d_CVI_shipment_rate_vec(ft, o_pcl, d_pcl, time_period,
//                                                                                                 up_map); //By reveiving county.
//                        rec_county_outcomes.clear();
//                        rec_county_outcomes.resize(allCounties.size(), 0);
//                        gsl_ran_multinomial(R, allCounties.size(), n_shipments,
//                                            d_lambdas.data(), rec_county_outcomes.data());
//
//                        //For every county that receives at least one shipment determine sending and rec premises.
//                        size_t n_processed_shipments_cc = 0; //County-county shipments.
//                        for(size_t d_county_idx=0; d_county_idx<allCounties.size(); ++d_county_idx)
//                        {
//                            if(rec_county_outcomes[d_county_idx] > 0)
//                            {
//                                size_t n_cc_shipments = rec_county_outcomes[d_county_idx]; //The number of shipments between this pair of counties.
//                                oprem_outcome.clear();
//                                oprem_outcome.resize(internal_prem_weights_o.size(), 0);
//                                gsl_ran_multinomial(R, internal_prem_weights_o.size(), n_cc_shipments,
//                                                    internal_prem_weights_o.data(), oprem_outcome.data());
//
//
//                                County* d_county = allCounties[d_county_idx];
//                                std::vector<Farm*> dcounty_premises = d_county->get_premises();
//                                const std::vector<double>& internal_prem_weights_d =
//                                    d_county->get_USAMMv3_dest_prem_weights_CVI(ft, d_pcl, time_period, up_map);
//                                size_t n_processed_shipments_pp = 0; //Prem-prem shipments.
//                                for(size_t oprem_idx=0; oprem_idx<oprem_outcome.size(); ++oprem_idx)
//                                {
//                                    if(oprem_outcome[oprem_idx] > 0)
//                                    {
//                                        Farm* oprem = ocounty_premises.at(oprem_idx);
//                                        size_t n_oprem_shipments = oprem_outcome[oprem_idx];
//                                        std::vector<unsigned int> dprem_outcome(internal_prem_weights_d.size(), 0);
//                                        gsl_ran_multinomial(R, internal_prem_weights_d.size(), n_oprem_shipments,
//                                                            internal_prem_weights_d.data(), dprem_outcome.data());
//
//                                        for(size_t dprem_idx=0; dprem_idx<dprem_outcome.size(); ++dprem_idx)
//                                        {
//                                            if(dprem_outcome[dprem_idx] > 0)
//                                            {
//                                                unsigned int n_pp_shipments = dprem_outcome[dprem_idx];
//                                                Farm* dprem = dcounty_premises.at(dprem_idx);
//                                                size_t shipment_volume;
//                                                double shipsize_nu = up.get_kNu(o_pcl, d_pcl, time_period);
//                                                double shipsize_mu = up.get_kMu(o_pcl, d_pcl, time_period);
//                                                for(unsigned int shipm_idx=0; shipm_idx<n_pp_shipments; ++shipm_idx)
//                                                {
//                                                                                    //Shipment between oprem and dprem.
//                                                    if(o_pcl == mkt_pcl and d_pcl == mkt_pcl)
//                                                    {
//                                                        //Both origin and receiver are markets, shipment size is a gamma-Poisson mixture RV.
//                                                        double mkt_shape = shipsize_nu; //Yes, mkt shipsize shape it's stored in the Nu map.
//                                                        double mkt_scale = shipsize_mu / mkt_shape; // mean / shape = scale
//                                                        double mkt_shipment_rate = draw_gamma(mkt_shape, mkt_scale);
//                                                        shipment_volume = draw_poisson(mkt_shipment_rate);
//                                                    }
//                                                    else
//                                                    {
//                                                        int prem_size_for_shipsize;
//                                                        double shipsize_alpha = shipsize_mu * shipsize_nu;
//                                                        double shipsize_beta = (1 - shipsize_mu) * shipsize_nu;
//                                                        double shipsize_p = draw_beta(shipsize_alpha, shipsize_beta); //The probability used in the binomial draw is a beta RV.
//                                                        if(o_pcl == mkt_pcl)
//                                                        {
//                                                            //Only sender is market, the shipment size is a beta-binomial RV where N is the receivers prem size.
//                                                            prem_size_for_shipsize = dprem->get_USAMMv3_binned_size();
//                                                        }
//                                                        else
//                                                        {
//                                                            //Neither premises is a market, shipment size is a beta-binomial RV where N is the senders prem size.
//                                                            prem_size_for_shipsize = oprem->get_USAMMv3_binned_size();
//                                                        }
//                                                        shipment_volume = draw_binom(prem_size_for_shipsize, shipsize_p);
//                                                    }
//                                                    Shipment* s = new Shipment{int(timestep), // timestep of shipment
//                                                                               day_of_year,
//                                                                               oprem,
//                                                                               dprem,
//                                                                               o_county->get_id(),
//                                                                               d_county->get_id(),
//                                                                               ft->get_species(),
//                                                                               shipment_volume,
//                                                                               true, //infectious
//                                                                               0, // ban (filled in filter_shipExposure if applicable)
//                                                                               time_period,
//                                                                               o_county->get_parent_state()->get_code(),
//                                                                               o_county->get_parent_state()->get_id(),
//                                                                               d_county->get_parent_state()->get_code(),
//                                                                               d_county->get_parent_state()->get_id(),
//                                                                               oprem->get_prem_class()->tag,
//                                                                               oprem->get_size_allSpecies(),
//                                                                               dprem->get_prem_class()->tag,
//                                                                               dprem->get_size_allSpecies()
//                                                                               };
//                                                    output.push_back(s);
//                                                }
//                                            }
//                                        }
//                                        n_processed_shipments_pp += n_oprem_shipments;
//                                        if(n_processed_shipments_cc == n_cc_shipments)
//                                        {
//                                            //All shipments from oprem have been taken care of,
//                                            //no use in checking the rest, they will all be 0.
//                                            break;
//                                        }
//                                    }
//                                }
//
//
//                                n_processed_shipments_cc += n_cc_shipments;
//                                if(n_processed_shipments_cc == n_shipments)
//                                {
//                                    //All shipments from o_county have been taken care of,
//                                    //no use in checking the rest, they will all be 0.
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

//void Shipment_manager::makeNetworkUSAMMv3(size_t timestep, size_t day_of_year,
//                                            std::string time_period, std::vector<Shipment*>& output,
//                                            std::vector<USAMMv3_parameters>& up_vec)
//{
//    std::vector<unsigned int> rec_county_outcomes; //To store the number of shipments to each county. Declared here so it can be reused.
//    output.clear();
//    rec_county_outcomes.reserve(allCounties.size());
//
//    //Determine the index of the market prem class (just for convenience).
//    Prem_class* mkt_pcl = nullptr;
//    std::vector<Prem_class*> pcl_vec = up_vec.begin()->getPremClasses();
//    for(Prem_class* pcl : pcl_vec)
//    {
//        if(pcl->tag == "Mkt")
//        {
//            mkt_pcl = pcl;
//            break;
//        }
//    }
//
//    for(auto& up : up_vec)
//    {
//        Farm_type* ft = up.getFarmType();
//        std::vector<Prem_class*> ft_prclasses = up.getPremClasses();
//        for(Prem_class* o_pcl : ft_prclasses)
//        {
//            for(Prem_class* d_pcl : ft_prclasses)
//            {
//                for(size_t o_county_idx=0; o_county_idx<allCounties.size(); ++o_county_idx)
//                {
//                    County* o_county = allCounties[o_county_idx];
//                    //This is the total rate of shipments from o_pcl to d_pcl leaving this county for any other county (or itself).
//                    double daily_ocounty_shipment_rate = o_county->get_o_CVI_shipment_rate(ft, o_pcl, d_pcl, time_period, up_vec);
//                    size_t n_shipments = gsl_ran_poisson(R, daily_ocounty_shipment_rate);
//
//                    if(n_shipments > 0)
//                    {
//                        //Determine destination counties.
//                        std::vector<double>& internal_prem_weights_o = o_county->get_USAMMv3_origin_prem_weights_CVI(ft, o_pcl, time_period, up_vec);
//                        std::vector<Farm*> ocounty_premises = o_county->get_premises();
//                        std::vector<unsigned int> oprem_outcome; //To store the number of shipments from different premises.
//                        oprem_outcome.reserve(internal_prem_weights_o.size());
//                        //These are the shipment rates from this county to each and every other county (incl. itself)
//                        //for shipments of the specific combination of prem classes. They are rates, but the probabilities
//                        //of the respective counties being receivers of the shipments are proportional to the rates. So
//                        //see them as unnormalized probabilities (weights). The multinomial function normalizes them internally.
//                        const std::vector<double>& d_lambdas = o_county->get_d_CVI_shipment_rate_vec(ft, o_pcl, d_pcl, time_period,
//                                                                                                 up_vec); //By reveiving county.
//                        rec_county_outcomes.clear();
//                        rec_county_outcomes.resize(allCounties.size(), 0);
//                        gsl_ran_multinomial(R, allCounties.size(), n_shipments,
//                                            d_lambdas.data(), rec_county_outcomes.data());
//
//                        //For every county that receives at least one shipment determine sending and rec premises.
//                        size_t n_processed_shipments_cc = 0; //County-county shipments.
//                        for(size_t d_county_idx=0; d_county_idx<allCounties.size(); ++d_county_idx)
//                        {
//                            if(rec_county_outcomes[d_county_idx] > 0)
//                            {
//                                size_t n_cc_shipments = rec_county_outcomes[d_county_idx]; //The number of shipments between this pair of counties.
//                                oprem_outcome.clear();
//                                oprem_outcome.resize(internal_prem_weights_o.size(), 0);
//                                gsl_ran_multinomial(R, internal_prem_weights_o.size(), n_cc_shipments,
//                                                    internal_prem_weights_o.data(), oprem_outcome.data());
//
//
//                                County* d_county = allCounties[d_county_idx];
//                                std::vector<Farm*> dcounty_premises = d_county->get_premises();
//                                const std::vector<double>& internal_prem_weights_d =
//                                    d_county->get_USAMMv3_dest_prem_weights_CVI(ft, d_pcl, time_period, up_vec);
//                                size_t n_processed_shipments_pp = 0; //Prem-prem shipments.
//                                for(size_t oprem_idx=0; oprem_idx<oprem_outcome.size(); ++oprem_idx)
//                                {
//                                    if(oprem_outcome[oprem_idx] > 0)
//                                    {
//                                        Farm* oprem = ocounty_premises.at(oprem_idx);
//                                        size_t n_oprem_shipments = oprem_outcome[oprem_idx];
//                                        std::vector<unsigned int> dprem_outcome(internal_prem_weights_d.size(), 0);
//
//                                        if(o_county == d_county and o_pcl == d_pcl)
//                                        {
//                                            std::vector<double> internal_prem_weights_d_copy = internal_prem_weights_d;
//                                            internal_prem_weights_d_copy[oprem->get_idx_in_county()] = 0.0; //Cannot ship to itself.
//                                            gsl_ran_multinomial(R, internal_prem_weights_d.size(), n_oprem_shipments,
//                                                                internal_prem_weights_d_copy.data(), dprem_outcome.data());
//                                        }
//                                        else
//                                        {
//                                            gsl_ran_multinomial(R, internal_prem_weights_d.size(), n_oprem_shipments,
//                                                                internal_prem_weights_d.data(), dprem_outcome.data());
//                                        }
//
//                                        for(size_t dprem_idx=0; dprem_idx<dprem_outcome.size(); ++dprem_idx)
//                                        {
//                                            if(dprem_outcome[dprem_idx] > 0)
//                                            {
//                                                unsigned int n_pp_shipments = dprem_outcome[dprem_idx];
//                                                Farm* dprem = dcounty_premises.at(dprem_idx);
//                                                for(unsigned int shipm_idx=0; shipm_idx<n_pp_shipments; ++shipm_idx)
//                                                {
//                                                    size_t shipment_volume;
//                                                    if(up.getUSAMMv3SpeciesCode() == "s")
//                                                    {
//                                                        shipment_volume = make_swine_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up);
//                                                    }
//                                                    else
//                                                    {
//                                                        shipment_volume = make_cattle_shipment_volume(oprem, dprem, o_pcl, d_pcl, time_period, up, mkt_pcl);
//                                                    }
//
//                                                    Shipment* s = new Shipment{int(timestep), // timestep of shipment
//                                                                               day_of_year,
//                                                                               oprem,
//                                                                               dprem,
//                                                                               o_county->get_id(),
//                                                                               d_county->get_id(),
//                                                                               ft->get_species(),
//                                                                               static_cast<int>(shipment_volume),
//                                                                               -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
//                                                                               0, // ban (filled in filter_shipExposure if applicable)
//                                                                               time_period,
//                                                                               o_county->get_parent_state()->get_code(),
//                                                                               o_county->get_parent_state()->get_id(),
//                                                                               d_county->get_parent_state()->get_code(),
//                                                                               d_county->get_parent_state()->get_id(),
//                                                                               oprem->get_prem_class()->tag,
//                                                                               oprem->get_size_allSpecies(),
//                                                                               dprem->get_prem_class()->tag,
//                                                                               dprem->get_size_allSpecies(),
//                                                                               nullptr,
//                                                                               "CVI"};
//                                                    output.push_back(s);
//                                                }
//                                            }
//                                        }
//                                        n_processed_shipments_pp += n_oprem_shipments;
//                                        if(n_processed_shipments_cc == n_cc_shipments)
//                                        {
//                                            //All shipments from oprem have been taken care of,
//                                            //no use in checking the rest, they will all be 0.
//                                            break;
//                                        }
//                                    }
//                                }
//
//
//                                n_processed_shipments_cc += n_cc_shipments;
//                                if(n_processed_shipments_cc == n_shipments)
//                                {
//                                    //All shipments from o_county have been taken care of,
//                                    //no use in checking the rest, they will all be 0.
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}

void Shipment_manager::makeNetworkUSAMMv2(std::vector<std::string> out_fnames, Grid_manager& G)
{
    std::vector<std::ofstream*> f_vec;
    for(std::string fname : out_fnames)
    {
        std::ofstream* f = new std::ofstream(fname + ".network");
        if(f->is_open())
        {
            //Write the header to output file.
            *f << "index\toCountyId\tdCountyId\tdayOfYear\tvolume\tunused\tperiod\toStateAbbr\toStateId\tdStateAbbr\tdStateId\toPremType\toPremSize\tdPremType\tdPremSize" <<
                  //"\tdistance\tunnormProb\tdState_inflow\tdState_n\tdCounty_n\tdCounty_weight" <<
                  std::endl;
            f_vec.push_back(f);
        }
        else
        {
            std::cout << "Failed to open network generation output file: " << fname
                      << ". Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }


    std::cout << "Generating shipments..." << std::endl;
    std::vector<Shipment*> shipments;
    shipments.reserve(1000000);
    size_t shipment_counter = 1;
    for(int day_i = 1; day_i < parameters->timesteps+1; day_i++)
    {
        int day_of_year = get_day_of_year(day_i, parameters->start_day);
        G.updateShippingParameters(day_i, day_of_year);
        std::string time_period = G.get_time_period();
        size_t rem_days = G.get_rem_days_of_period();
        std::vector<Farm*> dummy_vec = {};
        makeShipmentsUSAMMv2(day_i, rem_days, time_period, shipments, dummy_vec,
                             G.get_farm_types());
        for(Shipment* s : shipments)
        {
            Farm_type* ft = G.get_farm_type_by_name(s->species);
            size_t ft_index = ft->get_index();
            std::stringstream ss;
            ss << shipment_counter << "\t"
                        << s->origFIPS << "\t"
                        << s->destFIPS << "\t"
                        << s->day_of_year << "\t"
                        << -1 << "\t"
                        << -1 << "\t"
                        << time_period << "\t"
                        << s->origState_abbrev << "\t"
                        << s->origState_id << "\t"
                        << s->destState_abbrev << "\t"
                        << s->destState_id << "\t"
                        << s->originIndType << "\t"
                        << s->originSize << "\t"
                        << s->destIndType << "\t"
                        << s->destSize
//                        << 0.0 << "\t" //Distance
//                        << 0.0 << "\t" //Probability
//                        << 0.0 << "\t" //Dest state inflow
//                        << 0 << "\t" //Dest state n farms
//                        << 0 << "\t" //Dest county n farms
//                        << 0.0 << "\t" //Dest county inflow (weight).
//                        << s->origID << "\t" //Origin farm.
//                        << s->destID << "\t" //Destination farm.
                        << std::endl;
            *(f_vec.at(ft_index)) << ss.str();
            shipment_counter += 1;
            delete s;
        }
        shipments.clear();
    }
    for(size_t i = 0; i > f_vec.size(); i++)
    {
        f_vec[i]->close();
        delete f_vec[i];
    }

    std::vector<Farm_type*> ft_vec = G.get_farm_types();
    for(size_t i = 0; i < ft_vec.size(); i++)
    {
        std::ofstream f(out_fnames.at(i) + ".gen");
        if(f.is_open())
        {
            f << G.get_generation_string(ft_vec.at(i));
            f.close();
        }
    }
    std::cout << "...done." << std::endl;
}

void Shipment_manager::makeNetworkUSAMMv3(std::vector<std::string>out_fnames, Grid_manager& G)
{
    std::vector<std::ofstream*> f_vec;
    for(std::string fname : out_fnames)
    {
        std::ofstream* f = new std::ofstream(fname + ".network");
        if(f->is_open())
        {
            //Write the header to output file.
            *f << "index\toCountyId\tdCountyId\tdayOfYear\tvolume\tunused\tperiod\toStateAbbr\toStateId\tdStateAbbr\tdStateId\tpermitType\toPremType\toPremSize\tdPremType\tdPremSize\toPremId\tdPremId" <<
                  //"\tdistance\tunnormProb\tdState_inflow\tdState_n\tdCounty_n\tdCounty_weight" <<
                  std::endl;
            f_vec.push_back(f);
        }
        else
        {
            std::cout << "Failed to open network generation output file: " << fname
                      << ". Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }



    std::vector<Shipment*> shipments;
    shipments.reserve(1000000);
    size_t shipment_counter = 1;
    std::vector<Farm*> dummy_vec = {};
//    random_unique(G.get_allFarms_vector(), 300, dummy_vec);

//    for(int time_period_idx=0; time_period_idx<parameters->USAMM_temporal_n_timesteps.size(); ++time_period_idx)
//    {
//
//        std::string time_period = parameters->USAMM_temporal_order[time_period_idx];
//        G.updateShippingParameters(parameters->USAMM_temporal_start_points[time_period_idx],
//                                   parameters->USAMM_temporal_start_points[time_period_idx]);
//        makeShipmentsUSAMMv3(1, 1, time_period, time_period_idx,
//                             shipments, dummy_vec, G.get_usammv3_parameters_vec());
//        for(Shipment* s : shipments)
//        {
//            Farm_type* ft = G.get_farm_type_by_name(s->species);
//            size_t ft_index = ft->get_index();
//            *(f_vec.at(ft_index)) << shipment_counter << "\t"
//                        << s->origFIPS << "\t"
//                        << s->destFIPS << "\t"
//                        << s->day_of_year << "\t"
//                        << s->volume << "\t"
//                        << -1 << "\t"
//                        << time_period << "\t"
//                        << s->origState_abbrev << "\t"
//                        << s->origState_id << "\t"
//                        << s->destState_abbrev << "\t"
//                        << s->destState_id << "\t"
//                        << s->permitType << "\t"
//                        << s->originIndType << "\t"
//                        << s->originSize << "\t"
//                        << s->destIndType << "\t"
//                        << s->destSize << "\t"
//                        << s->oPrem->get_id() << "\t" //Origin farm.
//                        << s->dPrem->get_id() //Destination farm.
////                        << 0.0 << "\t" //Distance
////                        << 0.0 << "\t" //Probability
////                        << 0.0 << "\t" //Dest state inflow
////                        << 0 << "\t" //Dest state n farms
////                        << 0 << "\t" //Dest county n farms
////                        << 0.0 << "\t" //Dest county inflow (weight).
//
//                        << std::endl;
//            shipment_counter += 1;
//            delete s;
//        }
//        shipments.clear();
//    }

    for(int day_i = 1; day_i < parameters->timesteps+1; day_i++)
    {
        size_t day_of_year = get_day_of_year(day_i, parameters->start_day);
        G.updateShippingParameters(day_i, day_of_year);
        std::string time_period = G.get_time_period();
        int time_period_idx = G.get_temporal_index();

        makeShipmentsUSAMMv3(day_i, day_of_year, time_period, time_period_idx,
                             shipments, dummy_vec, G.get_usammv3_parameters_vec());
        for(Shipment* s : shipments)
        {
            Farm_type* ft = G.get_farm_type_by_name(s->species);
            size_t ft_index = ft->get_index();
            *(f_vec.at(ft_index)) << shipment_counter << "\t"
                        << s->origFIPS << "\t"
                        << s->destFIPS << "\t"
                        << s->day_of_year << "\t"
                        << s->volume << "\t"
                        << -1 << "\t"
                        << time_period << "\t"
                        << s->origState_abbrev << "\t"
                        << s->origState_id << "\t"
                        << s->destState_abbrev << "\t"
                        << s->destState_id << "\t"
                        << s->permitType << "\t"
                        << s->originIndType << "\t"
                        << s->originSize << "\t"
                        << s->destIndType << "\t"
                        << s->destSize << "\t"
                        << s->oPrem->get_id() << "\t" //Origin farm.
                        << s->dPrem->get_id() //Destination farm.
//                        << 0.0 << "\t" //Distance
//                        << 0.0 << "\t" //Probability
//                        << 0.0 << "\t" //Dest state inflow
//                        << 0 << "\t" //Dest state n farms
//                        << 0 << "\t" //Dest county n farms
//                        << 0.0 << "\t" //Dest county inflow (weight).

                        << std::endl;
            shipment_counter += 1;
            delete s;
        }
        shipments.clear();
    }

    for(size_t i = 0; i > f_vec.size(); i++)
    {
        f_vec[i]->close();
        delete f_vec[i];
    }

    std::vector<Farm_type*> ft_vec = G.get_farm_types();
    for(size_t i = 0; i < ft_vec.size(); i++)
    {
        std::ofstream f(out_fnames.at(i) + ".gen");
        if(f.is_open())
        {
            f << G.get_generation_string(ft_vec.at(i));
            f.close();
        }
    }
}

Shipment* Shipment_manager::generateShipmentUSAMMv2(Farm* origin_farm, size_t timestep, size_t day_of_year,
                                                    const std::string& time_period)
{
    //Generate a destination county until one is found that has the type required.
    Farm_type* origin_type = origin_farm->get_farm_type();
    County* origin_county = origin_farm->get_parent_county();
    County* dest_county = origin_county->get_shipment_destination_county(origin_type);

    while(dest_county->get_premises(origin_type).size() == 0)
    {
if(verbose>1){std::cout << "The destination county (" << dest_county->get_id() <<
                    ") does not contain any farms of type " << origin_type->get_species() <<
                    ", but has a probability of receiving shipments of that type. " <<
                    "This is likely due to the county having an inflow parameter != 0, " <<
                    "while simultaneously no premises of this type in the FLAPS data file " <<
                    "has this county as its parent. Choosing another..." << std::endl;}
        dest_county = origin_county->get_shipment_destination_county(origin_type);
    }

    //Pick one random element from the dest. countys vector of farms of correct type.
    Farm* destination_farm = dest_county->get_shipment_destination_premises(origin_type);
    int shipment_volume = -1;
    return new Shipment{static_cast<int>(timestep), // timestep of shipment
                        day_of_year,
                        origin_farm,
                        destination_farm,
                        origin_county->get_id(),
                        dest_county->get_id(),
                        origin_type->get_species(),
                        shipment_volume,
                        -1, //Number of infected animals on this shipment, unkown at this moment - determined in eval_exposure in and btb_eval_exposure in Status_manager
                        0, // ban (filled in filter_shipExposure if applicable)
                        time_period,
                        origin_county->get_parent_state()->get_code(),
                        origin_county->get_parent_state()->get_id(),
                        dest_county->get_parent_state()->get_code(),
                        dest_county->get_parent_state()->get_id(),
                        origin_farm->get_prem_class()->tag,
                        origin_farm->get_size_allSpecies(),
                        destination_farm->get_prem_class()->tag,
                        destination_farm->get_size_allSpecies(),
                        nullptr
                        };
}

size_t Shipment_manager::make_cattle_shipment_volume(Farm* o_prem, Farm* d_prem, Prem_class* o_pcl, Prem_class* d_pcl,
                                                     std::string time_period, USAMMv3_parameters& up, Prem_class* mkt_pcl)
{
    double shipsize_nu = up.get_kNu(o_pcl, d_pcl, time_period);
    double shipsize_mu = up.get_kMu(o_pcl, d_pcl, time_period);
    if(o_pcl == mkt_pcl and d_pcl == mkt_pcl)
    {
        //Both origin and receiver are markets, shipment size is a gamma-Poisson mixture RV.
        double mkt_shape = shipsize_nu; //Yes, mkt shipsize shape it's stored in the Nu map.
        double mkt_scale = shipsize_mu / mkt_shape; // mean / shape = scale
        double mkt_shipment_rate = draw_gamma(mkt_shape, mkt_scale);
        return draw_poisson(mkt_shipment_rate);
    }
    else
    {
        int prem_size_for_shipsize;
        double shipsize_alpha = shipsize_mu * shipsize_nu;
        double shipsize_beta = (1 - shipsize_mu) * shipsize_nu;
        double shipsize_p = draw_beta(shipsize_alpha, shipsize_beta); //The probability used in the binomial draw is a beta RV.
        if(o_pcl == mkt_pcl)
        {
            //Only sender is market, the shipment size is a beta-binomial RV where N is the receivers prem size.
            prem_size_for_shipsize = d_prem->get_USAMMv3_binned_size();
        }
        else
        {
            //Neither premises is a market, shipment size is a beta-binomial RV where N is the senders prem size.
            prem_size_for_shipsize = o_prem->get_USAMMv3_binned_size();
        }
        return draw_binom(prem_size_for_shipsize, shipsize_p);
    }
}

size_t Shipment_manager::make_swine_shipment_volume(Farm* o_prem, Farm* d_prem, Prem_class* o_pcl, Prem_class* d_pcl,
                                                    std::string time_period, USAMMv3_parameters& up)
{
    //All swine shipment sizes are gamma-Poisson with shape=swine_shipsize_shape and
    //mean = intercept * e^(osize*osize_scaling + dsize*dsize_scaling).

    double lognormed_osize = std::log(o_prem->get_USAMMv3_binned_size() / up.get_avg_prem_size(o_pcl));
    double lognormed_dsize = std::log(d_prem->get_USAMMv3_binned_size() / up.get_avg_prem_size(d_pcl));
    double oscaling;
    double dscaling;
    double shape;
    double intercept;

    if(o_prem->getCOMIdentifier() == nullptr) //CVI shipment (if one prem isn't a COM prem, the other won't be either).
    {
        oscaling = up.get_oShipSizeScalingCVI(o_pcl, d_pcl, time_period);
        dscaling = up.get_dShipSizeScalingCVI(o_pcl, d_pcl, time_period);
        shape = up.get_shipSizeShapeCVI(o_pcl, d_pcl, time_period);
        intercept = up.get_shipSizeInterceptCVI(o_pcl, d_pcl, time_period);

    }
    else
    {
        oscaling = up.get_oShipSizeScalingCOM(o_pcl, d_pcl, time_period);
        dscaling = up.get_dShipSizeScalingCOM(o_pcl, d_pcl, time_period);
        shape = up.get_shipSizeShapeCOM(o_pcl, d_pcl, time_period);
        intercept = up.get_shipSizeInterceptCOM(o_pcl, d_pcl, time_period);
    }

    double shipsize_mean = intercept * std::exp(lognormed_osize*oscaling + lognormed_dsize*dscaling);
    double rate = shape / shipsize_mean;
    double shipment_lambda_parameter = gsl_ran_gamma(R, shape, 1.0 / rate);
    return 1 + gsl_ran_poisson(R, shipment_lambda_parameter);
}

