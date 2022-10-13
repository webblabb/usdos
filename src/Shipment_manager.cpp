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
///	\param[in]	in_S			Pointer to Status_manager instance for this replicate
///	\param[in]	ffm				farm assignment method
///	\param[in]	speciesOnPrems	List of species
Shipment_manager::Shipment_manager(
	const std::vector<County*> in_FIPS_vec,
	const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpMap,
	Status_manager* in_S,
	const std::vector<std::string>& speciesOnPrems,
	const Parameters* p) :
        allCounties(in_FIPS_vec),
        fipsSpeciesMap(fipsSpMap),
        parameters(p),
        S(in_S),
        species(speciesOnPrems)
{
	verbose = verboseLevel;
	initialize();
}

Shipment_manager::Shipment_manager(
	const std::vector<County*> in_FIPS_vec,
	const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<Farm*> >>* fipsSpMap,
	const std::vector<std::string>& speciesOnPrems,
	const Parameters* p) :
        allCounties(in_FIPS_vec),
        fipsSpeciesMap(fipsSpMap),
        parameters(p),
        species(speciesOnPrems)
{
    S = nullptr;
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
        if(verbose>0){std::cout << "Shipment manager constructed: "<<allCounties.size()<<" counties with premises."
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

void Shipment_manager::makeShipmentsUSAMMv3(size_t timestep, size_t day_of_year, std::string time_period,
                                            int time_period_idx, std::vector<Shipment*>& output,
                                            std::vector<Farm*>& infFarms, std::vector<USAMMv3_parameters>& up_vec)
{
    int n_ships_generated = 0;

    //Select what farms will be involved in the generation of shipments.
    //Running with an empty infFarms is a signal to generating a full network of shipments,
    //so then we use all farms.
    //This function does have anything to do with btb, as btb is always simulated with USAMMv3
    std::vector<Farm*> affected_farms;
    affected_farms.reserve(1000000);
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

    //Determine the index of the market prem class (just for convenience/speed).
    Prem_class* mkt_pcl = nullptr;
    std::vector<Prem_class*> pcl_vec = up_vec.begin()->getPremClasses();
    for(Prem_class* pcl : pcl_vec)
    {
        if(pcl->tag == "Mkt")
        {
            mkt_pcl = pcl;
            break;
        }
    }


    //This function is only for making shipments during a disease simulation.
    for(Farm* oprem : affected_farms)
    {
        int n_farm_shipments = 0;
        double outgoing_shipment_rate = 0.0; //Keeps track of the sum of the rate of this prem sending to all other prems. For calculating the probability of sending to slaughter, used for the slaughter surveillance component.
        Prem_class* o_pcl = oprem->get_prem_class();
        std::vector<Farm_type*> o_ftypes;
        if(o_pcl == mkt_pcl)
        {
            //If it's a market, generate shipments for both beef and dairy.
            for(auto& up : up_vec)
            {
                o_ftypes.push_back(up.getFarmType());
            }
        }
        else
        {
            o_ftypes = { oprem->get_farm_type() };
        }

        for(Farm_type* o_fty : o_ftypes)
        {
            int o_fty_idx = o_fty->get_index();
            USAMMv3_parameters& up = up_vec[o_fty_idx];
            std::vector<Prem_class*> ft_pclasses = up.getPremClasses();
            County* o_county = oprem->get_parent_county();
            Shipment_kernel k(o_county->get_parent_state()->get_a(o_fty),
                              o_county->get_parent_state()->get_b(o_fty),
                              parameters->shipment_kernel, true);

            double normed_opremsize = oprem->get_USAMMv3_binned_size() / up.get_avg_prem_size(o_pcl);
            double ophi = up.get_phi_O(o_pcl, time_period);
            double oprem_size_w = std::pow(normed_opremsize, ophi);
            double outflow = o_county->get_parent_state()->get_outflow(o_fty);
            double ocov = o_county->get_county_ocov_weight(o_fty);

            for(County* d_county : allCounties)
            {
                for(Prem_class* d_pcl : d_county->get_prem_classes_by_type_idx(o_fty_idx)) //Only go through the pclasses that are present in the receiving county.
                {
                    double d_county_premsize_rep = d_county->get_rec_premsize_rep(o_fty, o_pcl, d_pcl, up_vec, time_period, time_period_idx); //Includes all prems in dcounty, their size, dphi and c.
                    if(o_county == d_county and o_pcl == d_pcl)
                    {
                        if(d_county_premsize_rep == 0)
                        {
                            continue; //I'm not entirely sure this will ever happen, but this check is not done very often so I'll leave it.
                        }
                        //Within county shipment needs to correct for the fact that a prem can't ship to itself.
                        d_county_premsize_rep -= std::pow(normed_opremsize, up.get_phi_D(o_pcl, time_period)) *
                                                 up.get_c(o_pcl, d_pcl, time_period);

                    }

                    double inflow = d_county->get_parent_state()->get_inflow(o_fty);
                    double dcov = d_county->get_county_dcov_weight(o_fty);
                    double k_val = k.kernel(o_county, d_county);
                    double d_county_receiving_rate = d_county_premsize_rep *
                                                     oprem_size_w *
                                                     outflow *
                                                     inflow *
                                                     ocov *
                                                     dcov *
                                                     k_val;
                    int n_county_shipments = draw_poisson(d_county_receiving_rate);
                    if(d_county_receiving_rate < 0.0)
                    {
                        std::cout << "Negative shipping rate." << std::endl;
                    }
                    outgoing_shipment_rate += d_county_receiving_rate;
                    if(n_county_shipments > 0)
                    {
                        int n_processed_shipments = 0;
                        std::vector<double> dcounty_internal_prem_weights = d_county->get_USAMMv3_dest_prem_weights(o_fty, d_pcl, time_period, up_vec);
                        std::vector<unsigned int> dprem_n_ships;
                        if(o_county == d_county and o_pcl == d_pcl)
                        {
                            //Remove the weight of the origin prem itself to prevent shipments to itself.
                            std::vector<double> dcounty_internal_prem_weights_copy = dcounty_internal_prem_weights;
                            dcounty_internal_prem_weights_copy[oprem->get_idx_in_county()] = 0.0;
                            draw_multinomial(n_county_shipments, dcounty_internal_prem_weights_copy, dprem_n_ships);
                        }
                        else
                        {
                            draw_multinomial(n_county_shipments, dcounty_internal_prem_weights, dprem_n_ships);
                        }

                        std::vector<Farm*> dprem_ptrs = d_county->get_premises();
                        for(size_t dprem_idx=0; dprem_idx<dprem_n_ships.size(); ++dprem_idx)
                        {
                            if(dprem_n_ships[dprem_idx] > 0)
                            {
                                Farm* dprem = dprem_ptrs[dprem_idx];
                                //Origin and destination premises determined, now determine shipment size.
                                size_t shipment_volume;
                                double shipsize_nu = up.get_kNu(o_pcl, d_pcl, time_period);
                                double shipsize_mu = up.get_kMu(o_pcl, d_pcl, time_period);
                                for(size_t ship_idx=0; ship_idx<dprem_n_ships[dprem_idx]; ++ship_idx)
                                {
                                    //Shipment between oprem and dprem.
                                    if(o_pcl == mkt_pcl and d_pcl == mkt_pcl)
                                    {
                                        //Both origin and receiver are markets, shipment size is a gamma-Poisson mixture RV.
                                        double mkt_shape = shipsize_nu; //Yes, mkt shipsize shape it's stored in the Nu map.
                                        double mkt_scale = shipsize_mu / mkt_shape; // mean / shape = scale
                                        double mkt_shipment_rate = draw_gamma(mkt_shape, mkt_scale);
                                        shipment_volume = 1 + draw_poisson(mkt_shipment_rate);
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
                                            prem_size_for_shipsize = dprem->get_USAMMv3_binned_size();
                                        }
                                        else
                                        {
                                            //Neither premises is a market, shipment size is a beta-binomial RV where N is the senders prem size.
                                            prem_size_for_shipsize = oprem->get_USAMMv3_binned_size();
                                        }
                                        shipment_volume = 1 + draw_binom(prem_size_for_shipsize - 1, shipsize_p);
                                    }

                                    Shipment* s = new Shipment{int(timestep), // timestep of shipment
                                                               day_of_year,
                                                               oprem,
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
                                                               oprem->get_prem_class()->tag,
                                                               oprem->get_size_allSpecies(),
                                                               dprem->get_prem_class()->tag,
                                                               dprem->get_size_allSpecies(),
                                                               nullptr};
                                    output.push_back(s);
                                    ++n_processed_shipments;
                                    ++n_ships_generated;
                                    ++n_farm_shipments;
                                    if(oprem == dprem)
                                    {
                                        std::cout << "Shipment origin == receiver." << std::endl;
                                    }
                                }
                            }
    //                        if(n_processed_shipments > 1000)
    //                        {
    //                            std::cout << d_county_premsize_rep << std::endl <<
    //                                         oprem_size_w << std::endl <<
    //                                                 outflow << std::endl <<
    //                                                 inflow << std::endl <<
    //                                                 ocov << std::endl <<
    //                                                 dcov << std::endl <<
    //                                                 k_val << std::endl <<
    //                                                 d_county_receiving_rate << std::endl;
    //                            std::cout << "fdewqfrwe" << std::endl;
    //                        }
                            if(n_processed_shipments >= n_county_shipments)
                            {
                                break;
                            }
                        }
                    }
                }
            }
        }
        oprem->set_latest_shipping_rate(outgoing_shipment_rate);
    } //End for farm.
//    std::cout << "Timestep " << timestep << ", " << n_ships_generated << " generated from " << infFarms.size() << " premises." << std::endl;
}

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
//                    double daily_ocounty_shipment_rate = o_county->get_o_shipment_rate(ft, o_pcl, d_pcl, time_period, up_map);
//                    size_t n_shipments = gsl_ran_poisson(R, daily_ocounty_shipment_rate);
//
//                    if(n_shipments > 0)
//                    {
//                        //Determine destination counties.
//                        std::vector<double>& internal_prem_weights_o = o_county->get_USAMMv3_origin_prem_weights(ft, o_pcl, time_period, up_map);
//                        std::vector<Farm*> ocounty_premises = o_county->get_premises();
//                        std::vector<unsigned int> oprem_outcome; //To store the number of shipments from different premises.
//                        oprem_outcome.reserve(internal_prem_weights_o.size());
//                        //These are the shipment rates from this county to each and every other county (incl. itself)
//                        //for shipments of the specific combination of prem classes. They are rates, but the probabilities
//                        //of the respective counties being receivers of the shipments are proportional to the rates. So
//                        //see them as unnormalized probabilities (weights). The multinomial function normalizes them internally.
//                        const std::vector<double>& d_lambdas = o_county->get_d_shipment_rate_vec(ft, o_pcl, d_pcl, time_period,
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
//                                    d_county->get_USAMMv3_dest_prem_weights(ft, d_pcl, time_period, up_map);
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

void Shipment_manager::makeNetworkUSAMMv3(size_t timestep, size_t day_of_year,
                                            std::string time_period, std::vector<Shipment*>& output,
                                            std::vector<USAMMv3_parameters>& up_vec)
{
    std::vector<unsigned int> rec_county_outcomes; //To store the number of shipments to each county. Declared here so it can be reused.
    output.clear();
    rec_county_outcomes.reserve(allCounties.size());

    //Determine the index of the market prem class (just for convenience/speed).
    Prem_class* mkt_pcl = nullptr;
    std::vector<Prem_class*> pcl_vec = up_vec.begin()->getPremClasses();
    for(Prem_class* pcl : pcl_vec)
    {
        if(pcl->tag == "Mkt")
        {
            mkt_pcl = pcl;
            break;
        }
    }

    for(auto& up : up_vec)
    {
        Farm_type* ft = up.getFarmType();
        std::vector<Prem_class*> ft_prclasses = up.getPremClasses();
        for(Prem_class* o_pcl : ft_prclasses)
        {
            for(Prem_class* d_pcl : ft_prclasses)
            {
                for(size_t o_county_idx=0; o_county_idx<allCounties.size(); ++o_county_idx)
                {
                    County* o_county = allCounties[o_county_idx];
                    //This is the total rate of shipments from o_pcl to d_pcl leaving this county for any other county (or itself).
                    double daily_ocounty_shipment_rate = o_county->get_o_shipment_rate(ft, o_pcl, d_pcl, time_period, up_vec);
                    size_t n_shipments = gsl_ran_poisson(R, daily_ocounty_shipment_rate);

                    if(n_shipments > 0)
                    {
                        //Determine destination counties.
                        std::vector<double>& internal_prem_weights_o = o_county->get_USAMMv3_origin_prem_weights(ft, o_pcl, time_period, up_vec);
                        std::vector<Farm*> ocounty_premises = o_county->get_premises();
                        std::vector<unsigned int> oprem_outcome; //To store the number of shipments from different premises.
                        oprem_outcome.reserve(internal_prem_weights_o.size());
                        //These are the shipment rates from this county to each and every other county (incl. itself)
                        //for shipments of the specific combination of prem classes. They are rates, but the probabilities
                        //of the respective counties being receivers of the shipments are proportional to the rates. So
                        //see them as unnormalized probabilities (weights). The multinomial function normalizes them internally.
                        const std::vector<double>& d_lambdas = o_county->get_d_shipment_rate_vec(ft, o_pcl, d_pcl, time_period,
                                                                                                 up_vec); //By reveiving county.
                        rec_county_outcomes.clear();
                        rec_county_outcomes.resize(allCounties.size(), 0);
                        gsl_ran_multinomial(R, allCounties.size(), n_shipments,
                                            d_lambdas.data(), rec_county_outcomes.data());

                        //For every county that receives at least one shipment determine sending and rec premises.
                        size_t n_processed_shipments_cc = 0; //County-county shipments.
                        for(size_t d_county_idx=0; d_county_idx<allCounties.size(); ++d_county_idx)
                        {
                            if(rec_county_outcomes[d_county_idx] > 0)
                            {
                                size_t n_cc_shipments = rec_county_outcomes[d_county_idx]; //The number of shipments between this pair of counties.
                                oprem_outcome.clear();
                                oprem_outcome.resize(internal_prem_weights_o.size(), 0);
                                gsl_ran_multinomial(R, internal_prem_weights_o.size(), n_cc_shipments,
                                                    internal_prem_weights_o.data(), oprem_outcome.data());


                                County* d_county = allCounties[d_county_idx];
                                std::vector<Farm*> dcounty_premises = d_county->get_premises();
                                const std::vector<double>& internal_prem_weights_d =
                                    d_county->get_USAMMv3_dest_prem_weights(ft, d_pcl, time_period, up_vec);
                                size_t n_processed_shipments_pp = 0; //Prem-prem shipments.
                                for(size_t oprem_idx=0; oprem_idx<oprem_outcome.size(); ++oprem_idx)
                                {
                                    if(oprem_outcome[oprem_idx] > 0)
                                    {
                                        Farm* oprem = ocounty_premises.at(oprem_idx);
                                        size_t n_oprem_shipments = oprem_outcome[oprem_idx];
                                        std::vector<unsigned int> dprem_outcome(internal_prem_weights_d.size(), 0);

                                        if(o_county == d_county and o_pcl == d_pcl)
                                        {
                                            std::vector<double> internal_prem_weights_d_copy = internal_prem_weights_d;
                                            internal_prem_weights_d_copy[oprem->get_idx_in_county()] = 0.0; //Cannot ship to itself.
                                            gsl_ran_multinomial(R, internal_prem_weights_d.size(), n_oprem_shipments,
                                                                internal_prem_weights_d_copy.data(), dprem_outcome.data());
                                        }
                                        else
                                        {
                                            gsl_ran_multinomial(R, internal_prem_weights_d.size(), n_oprem_shipments,
                                                                internal_prem_weights_d.data(), dprem_outcome.data());
                                        }

                                        for(size_t dprem_idx=0; dprem_idx<dprem_outcome.size(); ++dprem_idx)
                                        {
                                            if(dprem_outcome[dprem_idx] > 0)
                                            {
                                                unsigned int n_pp_shipments = dprem_outcome[dprem_idx];
                                                Farm* dprem = dcounty_premises.at(dprem_idx);
                                                size_t shipment_volume;
                                                double shipsize_nu = up.get_kNu(o_pcl, d_pcl, time_period);
                                                double shipsize_mu = up.get_kMu(o_pcl, d_pcl, time_period);
                                                for(unsigned int shipm_idx=0; shipm_idx<n_pp_shipments; ++shipm_idx)
                                                {
                                                                                    //Shipment between oprem and dprem.
                                                    if(o_pcl == mkt_pcl and d_pcl == mkt_pcl)
                                                    {
                                                        //Both origin and receiver are markets, shipment size is a gamma-Poisson mixture RV.
                                                        double mkt_shape = shipsize_nu; //Yes, mkt shipsize shape it's stored in the Nu map.
                                                        double mkt_scale = shipsize_mu / mkt_shape; // mean / shape = scale
                                                        double mkt_shipment_rate = draw_gamma(mkt_shape, mkt_scale);
                                                        shipment_volume = draw_poisson(mkt_shipment_rate);
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
                                                            prem_size_for_shipsize = dprem->get_USAMMv3_binned_size();
                                                        }
                                                        else
                                                        {
                                                            //Neither premises is a market, shipment size is a beta-binomial RV where N is the senders prem size.
                                                            prem_size_for_shipsize = oprem->get_USAMMv3_binned_size();
                                                        }
                                                        shipment_volume = draw_binom(prem_size_for_shipsize, shipsize_p);
                                                    }
                                                    Shipment* s = new Shipment{int(timestep), // timestep of shipment
                                                                               day_of_year,
                                                                               oprem,
                                                                               dprem,
                                                                               o_county->get_id(),
                                                                               d_county->get_id(),
                                                                               ft->get_species(),
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
                                                                               nullptr
                                                                               };
                                                    output.push_back(s);
                                                }
                                            }
                                        }
                                        n_processed_shipments_pp += n_oprem_shipments;
                                        if(n_processed_shipments_cc == n_cc_shipments)
                                        {
                                            //All shipments from oprem have been taken care of,
                                            //no use in checking the rest, they will all be 0.
                                            break;
                                        }
                                    }
                                }


                                n_processed_shipments_cc += n_cc_shipments;
                                if(n_processed_shipments_cc == n_shipments)
                                {
                                    //All shipments from o_county have been taken care of,
                                    //no use in checking the rest, they will all be 0.
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
            *f << "index\toCountyId\tdCountyId\tdayOfYear\tvolume\tunused\tperiod\toStateAbbr\toStateId\tdStateAbbr\tdStateId\toPremType\toPremSize\tdPremType\tdPremSize\toPremId\tdPremId" <<
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
        size_t day_of_year = get_day_of_year(day_i, parameters->start_day);
        G.updateShippingParameters(day_i, day_of_year);
        std::string time_period = G.get_time_period();
        std::vector<Farm*> dummy_vec = {};
        makeNetworkUSAMMv3(day_i, day_of_year, time_period, shipments,
                           G.get_usammv3_parameters_vec());
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
    std::cout << "...done." << std::endl;
}

Farm* Shipment_manager::largestStatus(std::vector<Farm*>& premVec, std::string& status)
{
    if(S == nullptr)
    {
        std::cout << "This shipment manager was created without status manager. "
                  << "Use the correct constructor." << std::endl << "Exiting..."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
	bool found = 0;
	auto i = premVec.back(); // start at end of sorted vector (largest prem) and work backwards
	while (found==0){
		if(S->getAny_diseaseStatus(i).compare(status)==0){ // if prem has this status, stop and return this prem
			found = 1;
		} else if ( i>premVec.front() ){i--; // keep moving backwards
		} else if ( i==premVec.front() ){ i = premVec.back(); found = 1;} // if no prems have this status, return largest
	}
	return i;
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
