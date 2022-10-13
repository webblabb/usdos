#include <stdio.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include "Farm.h"
#include "County.h"
#include "shared_functions.h"
#include "Shipment_manager.h"
#include "State.h"
#include "File_manager.h"

size_t Farm::current_quarter_idx = 0;

Farm::Farm(int in_id, double in_x, double in_y, std::string in_fips,
           const Parameters* parameters)
	:
	id(in_id),
	cellID(-1),
	x_coordinate(in_x),
	y_coordinate(in_y),
	position(in_x, in_y),
	fips(in_fips),
	parameters(parameters),
	btb_p(&parameters->btbWithinHerdParams),
	neighborRadiusCalculated(0)
{
    speciesCountsMapQuarterly.resize(4, std::unordered_map<std::string, int>());
    sus.resize(4, 0.0);
    inf.resize(4, 0.0);
    inf_max = 0.0;
    sus_max = 0.0;
}

Farm::~Farm()
{
}

int Farm::get_prem_class_idx() const
{
     return prem_class->idx;
}

int Farm::get_USAMMv3_unbinned_size() const
{
    return this->get_size_specific_quarter(this->farm_type, USAMMv3_unbinned_size_idx);
}

double Farm::get_normalized_oweight(Farm_type* ft) const
{
    return normalized_oweight.at(ft);
}

double Farm::get_unnormalized_oweight(Farm_type* ft) const
{
    return oweight.at(ft);
}

double Farm::get_unnormalized_dweight(Farm_type* ft) const
{
    return dweight.at(ft);
}

double Farm::get_shipment_rate()
{
    if(latest_shipping_rate < 0)
    {
        std::cout << "Attempting to get shipment rate from a premises "
                  << "without first updating it's shipping rate. Exiting." << std::endl;
        exit(EXIT_FAILURE);
    }

    return latest_shipping_rate;
}

const std::unordered_map<std::string, int>& Farm::get_spCountsCurrentQuarter()
{
    return speciesCountsMapQuarterly.at(Farm::current_quarter_idx);
}

const std::unordered_map<std::string, int>& Farm::get_spCountsSpecificQuarter(int quarter_idx)
{
    return speciesCountsMapQuarterly.at(quarter_idx);
}

int Farm::get_size_current_quarter(const Farm_type* ft) const
{
    return speciesCountsQuarterly[ft->get_index()][Farm::current_quarter_idx];
}

int Farm::get_size_specific_quarter(const Farm_type* ft, int quarter_idx) const
{
    return speciesCountsQuarterly[ft->get_index()][quarter_idx];
}

int Farm::get_size_allSpecies() const
{
	int num = 0;
	for (auto& i : speciesCountsQuarterly){
		num += i.at(Farm::current_quarter_idx);
	}
  return num;
}

void Farm::set_xy(double x, double y)
{
    x_coordinate = x;
    y_coordinate = y;
}

void Farm::set_cellID(const int in_cellID)
{
	cellID = in_cellID;
}

void Farm::set_farm_type(Farm_type* in_type)
{
    farm_type = in_type;
    farm_type_str = in_type->get_species();
}

void Farm::set_prem_class(Prem_class* in_class)
{
    prem_class = in_class;
    if(prem_class->tag == "Frm")
    {
        is_frm = true;
    }
    else if(prem_class->tag == "Fdl")
    {
        is_fdl = true;
    }
    else if(prem_class->tag == "Mkt")
    {
        is_mkt = true;
    }
}

void Farm::set_quarterlySpeciesCounts(int ft_idx, std::string ft_name, std::vector<int> sp_count)
{
    if(int(speciesCountsQuarterly.size()) < ft_idx+1)
    {
        speciesCountsQuarterly.emplace_back(sp_count);
    }
    else
    {
        speciesCountsQuarterly[ft_idx] = sp_count;
    }
    for(size_t i=0; i<sp_count.size(); ++i)
    {
        speciesCountsMapQuarterly[i][ft_name] = sp_count[i];
    }
}

void Farm::set_sus(const double in_sus, unsigned int quarter_idx)
{
	sus.at(quarter_idx) = in_sus;
	if(in_sus > sus_max)
    {
        sus_max = in_sus;
    }
}

void Farm::set_sus_all_quarters(const double in_sus)
{
    for(unsigned int quarter_idx=0; quarter_idx<4; ++quarter_idx)
    {
        this->set_sus(in_sus, quarter_idx);
    }
}

void Farm::set_inf(const double in_inf, unsigned int quarter_idx)
{
	inf.at(quarter_idx) = in_inf;
}

void Farm::set_inf_all_quarters(const double in_inf)
{
    for(unsigned int quarter_idx=0; quarter_idx<4; ++quarter_idx)
    {
        this->set_inf(in_inf, quarter_idx);
    }
}

void Farm::set_inf_max(const double in_inf_max)
{
    inf_max = in_inf_max;
}

void Farm::set_normalized_oweight(const double in_normalized_oweight, Farm_type* ft)
{
    normalized_oweight[ft] = in_normalized_oweight;
}

void Farm::set_normalized_dweight(const double in_normalized_dweight, Farm_type* ft)
{
    normalized_dweight[ft] = in_normalized_dweight;
}

void Farm::set_unnormalized_oweight(const double in_oweight, Farm_type* ft)
{
    oweight[ft] = in_oweight;
}

void Farm::set_unnormalized_dweight(const double in_dweight, Farm_type* ft)
{
    dweight[ft] = in_dweight;
}

void Farm::set_parent_county(County* in_county)
{
  parent_county = in_county;
}

/// Sets distances-neighbors map via swap
void Farm::set_distancesNeighbors(std::multimap<double, Farm*>& dn)
{
	distancesNeighbors.swap(dn);
}

State* Farm::get_parent_state() const
{
    return parent_county->get_parent_state();
}

void Farm::set_neighborRadiusCalculated(const double radius)
{
	neighborRadiusCalculated = radius;
}


Farm_type::Farm_type(std::string herd, std::vector<std::string> in_species, unsigned int index) :
    index(index),
    herd(herd)
{
    std::stringstream ss;
    for(size_t i = 0; i < herd.size(); i++)
    {
        if(herd[i] != '0')
        {
           ss << in_species[i] << ", ";
        }
    }
    species = ss.str();
    species.pop_back(); // get rid of last space
    species.pop_back(); // get rid of last comma
}

Farm_type::~Farm_type()
{
}

Prem_status::Prem_status(Farm* f)
	:
	Farm(*f), // call copy constructor of Farm to fill in all other members
	fileStatus("notDangerousContact"),
	diseaseStatus("sus"),
	isVaccinated(false)
{
    unvaccinate();
	probPreventExposure.emplace_back(0.0);
	probPreventTransmission.emplace_back(0.0);

	btb_infection_classes.resize(6,0);
	btb_infection_classes[0] = get_size_current_quarter(farm_type);
	verify_premsize();
	btb_n_infectious = 0;

}

Prem_status::~Prem_status()
{
}

int Prem_status::get_start(std::string s) const
{
	int output = -1;
	if (start.count(s) > 0){
		output = start.at(s);
	}
	return output;
}

int Prem_status::get_end(std::string s) const
{
	int output = -1;
	if (end.count(s) > 0){
		output = end.at(s);
	}
	return output;
}

int Prem_status::get_testResult(std::string d_type) const
{
	int output = -1;
	if (testResult.count(d_type) > 0){
		output = testResult.at(d_type);
	}
	return output;
}

/// Returns the probability that exposure of a premises is prevented, based on the
/// effectiveness of currently effective control types. If multiple control types are
/// currently effective, the maximum effectiveness value is returned (no additive
/// effects are assumed).
double Prem_status::get_probPreventExposure() const
{
	auto ppe = std::max_element(probPreventExposure.begin(), probPreventExposure.end());
	double output = *ppe;
	if (output > 1){output = 1;
	} else if (output < 0){output = 0;}
	return output;
}

void Prem_status::rem_probPreventExposure(double eff)
{
	auto it = std::find(probPreventExposure.begin(), probPreventExposure.end(), eff);
	if (it != probPreventExposure.end()) {
		// swap the one to be removed with the last element
		std::vector<double>::iterator last = std::prev(probPreventExposure.end());
		std::iter_swap(it, last);
		// and remove the item at the back of the container
		probPreventExposure.pop_back();
	}
}
/// Returns the probability that transmission from a premises is prevented, based on the
/// effectiveness of currently effective control types. If multiple control types are
/// currently effective, the maximum effectiveness value is returned (no additive
/// effects are assumed).
double Prem_status::get_probPreventTransmission() const
{
	auto ppe = std::max_element(probPreventTransmission.begin(), probPreventTransmission.end());
	double output = *ppe;
	if (output > 1){output = 1;
	} else if (output < 0){output = 0;}
	return output;
}

void Prem_status::rem_probPreventTransmission(double eff)
{
	auto it = std::find(probPreventTransmission.begin(), probPreventTransmission.end(), eff);
	if (it != probPreventTransmission.end()) {
		// swap the one to be removed with the last element
		std::vector<double>::iterator last = std::prev(probPreventTransmission.end());
		std::iter_swap(it, last);
		// and remove the item at the back of the container
		probPreventTransmission.pop_back();
	}
}

void Prem_status::vaccinate(double efficacy)
{
    if(efficacy > 1.0 or efficacy < 0.0)
    {
        std::cout << "Vaccine efficacy must be between 0 and 1." << std::endl;
        exit(EXIT_FAILURE);
    }
    //Assumes that all commodities are vaccinated
    for(auto& sp_count_pair : this->get_spCountsCurrentQuarter())
    {
        if(sp_count_pair.second > 0)
        {
            int n_unvaxed = draw_binom(sp_count_pair.second, 1.0 - efficacy);
            currentUnvaccinatedPrevalence[sp_count_pair.first] = double(n_unvaxed) / double(sp_count_pair.second);
        }
    }
    isVaccinated = true;
}

void Prem_status::unvaccinate()
{
    for(const auto& sp_n_pair : this->get_spCountsCurrentQuarter())
    {
        currentUnvaccinatedPrevalence[sp_n_pair.first] = 1.0;
    }
    isVaccinated = false;
}

///Returns the number of unvaccinated animals, if the herd is vaccinated, it's
///the part of the herd for which the vaccine was inefficient; if the herd isn't
///vaccinated it's all animals.
std::unordered_map<std::string, int> Prem_status::get_currentSizeUnvaccinated()
{
    if(isVaccinated)
    {
        std::unordered_map<std::string, int> unvaccinated_n = this->get_spCountsCurrentQuarter();
        for(auto& sp_num_pair : unvaccinated_n)
        {
            sp_num_pair.second = sp_num_pair.second * currentUnvaccinatedPrevalence[sp_num_pair.first];
        }
        return unvaccinated_n;
    }
    return this->get_spCountsCurrentQuarter();
}


/// Stores info for determining if potential DCs will be dangerousContacts
void Prem_status::add_potentialDCInfo(Farm* pdc, const std::vector<bool>& dcInfo)
{
	if (potentialDCs.count(pdc)==0){
		potentialDCs[pdc] = dcInfo;
		if(verboseLevel>1){std::cout<<"F in Prem_status add_potentialDCInfo first if statement"<<std::endl;}

	} else { // in case DC relationship already established between this pair of farms,
	         // isDC status of 'true' should override any 'false' evaluations
		if(verboseLevel>1){std::cout<<"F in Prem_status add_potentialDCInfo first if statement else"<<std::endl;}

if(verboseLevel>1){std::cout<<"DC relationship already established, new values: ";}
        std::vector<bool> existingPDC = potentialDCs.at(pdc);
        std::vector<DiseaseStatus> statuses;
        if(parameters->infectionType == InfectionType::BTB)
        {
            statuses = {DiseaseStatus::SUS, DiseaseStatus::BTB};
        }
        else
        {
            statuses = {DiseaseStatus::SUS, DiseaseStatus::EXP, DiseaseStatus::INF};
        }

        for(DiseaseStatus status : statuses)
        {
            int status_idx = static_cast<int>(status);
            if(dcInfo.at(status_idx))
            {
                potentialDCs.at(pdc).at(status_idx) = true;
                if(verboseLevel>1){std::cout<<"F in Prem_status add_potentialDCInfo for loop"<<std::endl;}
            }
            if(verboseLevel>1){std::cout<<parameters->DiseaseStatusToString.at(status)<<"="<<potentialDCs.at(pdc).at(status_idx)<<", ";}
        }
if(verboseLevel>1){std::cout<<std::endl;}
	}
}
/// Records time, source premises of exposure, route (local vs shipping), and whether or not
/// exposure was prevented.
void Prem_status::add_exposureSource(int t, Farm* source, int route, std::string block)
{
	// check that vectors are all the same size to ensure indices match up
	if (expTime.size() == expSource.size() &&
			expSource.size() == expRoute.size() &&
			expRoute.size() == expBlocked.size()){
		expTime.emplace_back(t);
		expSource.emplace_back(source);
		expRoute.emplace_back(route);
		expBlocked.emplace_back(block);
	} else {
		std::cout<<"ERROR: In Prem_status::add_exposureSource, exposureSource vectors are different sizes. Exiting...";
		exit(EXIT_FAILURE);
	}
}

/// Checks if premises f is the source of exposure. Used for Dangerous Contacts evaluation.
bool Prem_status::is_exposureSource(Farm* f)
{
	// expBlocked is generally very short - mostly length=1
	bool answer = 0;
	for (unsigned int s=0; s<expSource.size(); ++s){
		if (expBlocked.at(s).compare("none")==0){ // if exposure was not prevented
				if (expSource.at(s) == f){ // check if farm matches f
					answer = 1;
				}
		}
	}
	return answer;
}

/// Returns true if premises was not reported at time of addition to diagnostic waitlist
bool Prem_status::was_dcWhenDiagnosticWaitlisted(const std::string d_type){
	return (statusWhenDiagnosticWaitlisted.at(d_type) == "dangerousContact" || statusWhenDiagnosticWaitlisted.at(d_type) == "exposed" || statusWhenDiagnosticWaitlisted.at(d_type) == "suspected");
}

/// Returns true if premises was not reported at time of addition to waitlist
bool Prem_status::was_dcWhenWaitlisted(const std::string c_type){
	return (statusWhenWaitlisted.at(c_type) == "dangerousContact" || statusWhenWaitlisted.at(c_type) == "exposed");
}

/// Returns the time of infection for a farm
int  Prem_status::when_infected()
{
    int exposed_at = 999999;
    // must find the first time a farm was infected and the exposure was not prevented
    // expBlocked is generally very short - mostly length=1
    for (unsigned int s=0; s<expSource.size(); ++s){
        // If exposure was not prevented. Should happen at least once since this prem is infectious.
        // And if the time of the occurrece was before current value of exposed_at.
        if (expBlocked.at(s).compare("none")==0 and
            expTime.at(s) < exposed_at){
                exposed_at = expTime.at(s);
        }
    }

    // If exposed_at is unchanged since declaration, the exposure time has not been found.
    // source exposure time ==1
    if(exposed_at==999999){
        std::cout<<"ERROR: In Prem_status::when_infected, time of exposure not found. Exiting...";
        exit(EXIT_FAILURE);
    }
    return exposed_at;
}

///calculates current infectiousness for farm and current time t using
///the entire population size.
double Prem_status::get_inf_partial_as_unvaccinated(int t, int quarter_idx, const Parameters* p,
                                                    const std::unordered_map<std::string, std::vector<double>>& normInf_map)
{
    return get_inf_partial(t, quarter_idx, p, normInf_map, this->get_spCountsCurrentQuarter());
}

///calculates current infectiousness for farm and current time t assuming
///that the farm has been vaccinated and using the part of population size for which
///vaccination was unsuccessful.
double Prem_status::get_inf_partial_as_vaccinated(int t, int quarter_idx, const Parameters* p,
                                                  const std::unordered_map<std::string, std::vector<double>>& normInf_map)
{
    return get_inf_partial(t, quarter_idx, p, normInf_map, this->get_currentSizeUnvaccinated());
}

///calculates current infectiousness for farm and current time t
double Prem_status::get_inf_partial(int t, int quarter_idx, const Parameters* p,
                                    const std::unordered_map<std::string, std::vector<double>>& normInf_map,
                                    const std::unordered_map<std::string, int>& sp_counts)
{
    //If it's already been calculated for this timestep, return that.
    if(t == saved_evaluated_PTF_FMD_timestep)
    {
        return saved_evaluated_PTF_FMD;
    }

    //Otherwise, calculate and save for this timestep.
    //find time of infection
    int t0 = this->when_infected();
    int time_since_infection =  t - t0;
    double current_n_inf_PT_FMD = 0.0;
    double premInf = 0.0;
    //The following is the scaling factor for the PT-parameters.
    //For parameters r0=0.05, r1=0.006, g=0.44, ts0=4 and infectious
    //period = 30 days it should be set to 6.272089. It's independent of prem size.
    double factor = p->partialParams.at(4);
    double latency_FMD = std::get<0>(p->latencyParams);
    for (auto& sp_N_pair : sp_counts)
    {
        double n_infectious = get_effective_n_infectious_animals(time_since_infection, sp_N_pair.second, p);
        current_n_inf_PT_FMD += n_infectious;
        double normInf_val=normInf_map.at(sp_N_pair.first).at(quarter_idx);
        double spInf = normInf_val*std::pow(factor*n_infectious, p->infExponents.at(sp_N_pair.first)); // transmissibility value for this species/type
        premInf += spInf; // add this species to the total for this premises
    }

    //If this premises no longer have at least one effectively infectious animal, force premature ending of it's status as infectious
    //by changing the end time of the current status to the next timestep. This will cause it to become removed at the beginning of
    //the next time step. Also, in order to make sure that small farms don't get their infectious period ended in the beginning
    //(when the PT function will evaluate to a small value since it's still ramping up), make sure that at least a week of
    //infectiousness has passed before forcing removal.
    if(current_n_inf_PT_FMD < 1.0 and time_since_infection > 6+latency_FMD)
    {
        this->set_end("inf", t+1);
    }
    saved_evaluated_PTF_FMD = premInf;
    saved_evaluated_PTF_FMD_timestep = t;
    return saved_evaluated_PTF_FMD;
}

double Prem_status::get_max_inf_partial(const Parameters* p, const std::unordered_map<std::string, double>& normInf_map,
                                        const std::unordered_map<std::string, int>& sp_counts)
{
    double max_prem_inf = 0.0;
    double latency_FMD = std::get<0>(p->latencyParams);
    double inf_period_FMD = std::get<0>(p->infectiousParams);
    for(int time_since_infection = int(latency_FMD);
        time_since_infection < inf_period_FMD + latency_FMD;
        ++time_since_infection)
    {
        double premInf = 0.0;
        //The following is the scaling factor for the PT-parameters.
        //For parameters r0=0.05, r1=0.006, g=0.44, ts0=4 and infectious
        //period = 30 days it should be set to 6.272089. It's independent of prem size.
        double factor = p->partialParams.at(4);

        for (auto& sp_N_pair : sp_counts)
        {
            double n_infectious = get_effective_n_infectious_animals(time_since_infection, sp_N_pair.second, p);
            double normInf_val=normInf_map.at(sp_N_pair.first);
            double spInf = normInf_val*std::pow(factor*n_infectious, p->infExponents.at(sp_N_pair.first)); // transmissibility value for this species/type
            premInf += spInf; // add this species to the total for this premises
        }
        if(premInf > max_prem_inf)
        {
            max_prem_inf = premInf;
        }
    }
    return max_prem_inf;
}


///calculates current infectiousness for farm and current time t
double Prem_status::get_effective_n_infectious_animals(int time_since_infection,
                                                       int n_animals, const Parameters* p)
{
    ////fixed parameters
    double r0 = p->partialParams.at(0);
    double r1 = p->partialParams.at(1);
    double g = p->partialParams.at(2);
    double tS0 = p->partialParams.at(3);


    double I1 = 0.0;
    double I2 = 0.0;

    //tsigma is the same as the latency duration
    double tsigma = std::get<0>(p->latencyParams);

    //duration of time infected
    //must be converted to double for below calculation
    double tInf = double(time_since_infection) + 1.0; //This is used locally within this function and needs the +1 to work correctly with the latency between exp and inf as handled insude the status_manager

    double NInf=0.0; //Effective number of infectious animals, the result of the PT function.
    if(tInf>tsigma)
    {
        double N = double(n_animals);
        if(N == 0)
        {
            return 0.0;
        }

        ////calculate growth rate
        double r = N*(r0*(tInf-tsigma-1) + r1*std::pow(tInf-tsigma,2));

        ////evaluate growth curve
        if(tInf<(tsigma+tS0))
        {
            NInf += I1 +
                    I2*(tInf-tsigma) / tS0 +
                    (r + g*r*(tS0-(tInf-tsigma))) / (tS0*g*g) -
                    std::exp(-g*(tInf-tsigma)) * (r+g*r*tS0) / (tS0*g*g);
        }
        else if(tInf>=(tsigma+tS0))
        {
            NInf += I1 +
                    I2 +
                    std::exp(-g*(tInf-tsigma)) * (r*std::exp(g*tS0)-r-g*r*tS0) / (tS0*g*g);
        }

        if(N!=0&&NInf>N){
           std::cout<<"ERROR: In Grid_checker::get_inf_partial, the number of infectious animals is greater than the number of animals on the farm. Exiting...";
            exit(EXIT_FAILURE);
        }
        if(NInf < 0)
        {
            std::cout<<"ERROR: In Grid_checker::get_inf_partial, the number of infectious animals is less than zero. Exiting...";
            exit(EXIT_FAILURE);
        }
    }
    return NInf;
}

int Prem_status::btb_seed()
{
    //Assign 10% of the herd to the first exposed class.
    //This is just arbitrary at the moment.
    int n_to_seed = int(std::ceil(0.1 * this->get_spCountsCurrentQuarter().at(farm_type->get_species())));
    btb_infection_classes[0] -= n_to_seed;
    btb_infection_classes[1] += n_to_seed;
    return n_to_seed;
}

//void Prem_status::btb_initialize_false_pos()
//{
//    btb_true_neg = draw_binom(btb_infection_classes[0], btb_p->test_specificity * btb_p->specificity_scaling_factor);
//    btb_false_pos = btb_infection_classes[0] - btb_true_neg;
//}

bool Prem_status::btb_update_size_change_quarter()
{
    //Updates the number of animals in each infection class following a change
    //in prem size (e.g. due to different quarterly prem sizes).
    //If there are more animals now, add them all to the susceptible class,
    //if there are fewer, remove from all classes equaly.
    double current_size = double( std::accumulate(btb_infection_classes.begin(),
                                                  btb_infection_classes.end(), 0) );
    int prem_size = this->get_size_current_quarter(this->farm_type); //This changes based on quarter.

    int size_change = int(prem_size - current_size);

    if(size_change > 0) //Positive = add animals to sus.
    {
        btb_infection_classes[0] += size_change;
    }
    else if(size_change < 0) //Negative = remove equally from all.
    {
//        std::vector<double> removal_p(btb_infection_classes.size(), 1.0 / double(btb_infection_classes.size())); //Equal prob, for each class.
        std::vector<unsigned int> remove_n(btb_infection_classes.size(), 0);
        draw_multivariate_hypergeometric(std::abs(size_change), btb_infection_classes, remove_n);
        for(size_t i=0; i<btb_infection_classes.size(); ++i)
        {
            btb_infection_classes[i] -= remove_n[i];
        }
    }

    //If size change == 0 no action neeed.

    btb_n_infectious = btb_infection_classes[4] + btb_infection_classes[5];
    verify_premsize();
    if(int(btb_infection_classes[0]) == prem_size)
    {
        return false; //All animals are in the susceptible class, so return false to indicate no longer infected.
    }
    return true;
}

double Prem_status::btb_get_prevalence_infectious()
{
    verify_premsize();
    //Used as a measure of the premises infectiousness in the btb local spread component.
    size_t prem_size = this->get_spCountsCurrentQuarter().at(farm_type->get_species());
    double prevalence_inf = double(btb_n_infectious) / double(prem_size);
    return prevalence_inf;
}

int Prem_status::btb_get_N_infected()
{
    return std::accumulate(btb_infection_classes.begin()+1, btb_infection_classes.end(), 0);
//    return this->get_spCountsCurrentQuarter().at(farm_type->get_species()) - btb_infection_classes[0];
}

void Prem_status::btb_expose_susceptibles(unsigned int n_to_expose)
{
    verify_premsize();
    unsigned int max_to_expose = std::min(n_to_expose, btb_infection_classes[0]); //Can't expose more than whats in the susc class.
    btb_infection_classes[0] -= max_to_expose;
    btb_infection_classes[1] += max_to_expose;
    btb_n_new_external_infections += max_to_expose;
    verify_premsize();
}

void Prem_status::btb_add_to_inf_classes(const std::vector<unsigned int>& n_in_btb_classes)
{
    /*
    //Shipment brings in various numbers of animals to the different classes.
    //Since we don't model shipments from non-infected sources, we don't care about
    //susceptible animals on the shipment, but add only the animals in the other
    //classes and let them change the prevalence on the premises.
    unsigned int prem_size = speciesCounts[farm_type->get_species()]; //"True" size.
    verify_premsize();

    //Grow the classes with the additional infected animals.
    unsigned int new_size = btb_infection_classes[0];
    for(size_t i=1; i<btb_infection_classes.size(); ++i) //for all non-susc classes.
    {
        btb_infection_classes[i] += n_in_btb_classes[i];
        new_size += btb_infection_classes[i]; //Also keep track of the new total size.
    }

    //Now remove animals to get back to the original size, but preserving the prevalence.
    unsigned int n_to_remove = new_size - prem_size;
    if(n_to_remove > 0)
    {
        //Select animals to remove by making a sample from a multiv. hypergeometric distr.
        std::vector<unsigned int> removals;
        draw_multivariate_hypergeometric(n_to_remove, btb_infection_classes, removals);

        for(size_t i=0; i<btb_infection_classes.size(); ++i)
        {
            btb_infection_classes[i] -= removals[i];
        }
    }
    btb_n_infectious = btb_infection_classes[4] + btb_infection_classes[5];
    verify_premsize();
    */

    /////////////////////////////////////////////////////

    //Older version of this function that just ignores that prevalence can drop
    //significantly if a large shipment with low prevalence comes in and replaces
    //the "old" animals.

    //Shipment brings in various numbers of animals to the different classes.
    //Add them and then shrink them so that total vol remains the same, but with
    //new prevalences.
    unsigned int prem_size = this->get_spCountsCurrentQuarter().at(farm_type->get_species()); //"True" size.
    unsigned int current_size = std::accumulate(btb_infection_classes.begin(),
                                                btb_infection_classes.end(), 0);
    verify_premsize();

    //Grow the classes with the additional animals.
    for(size_t i=0; i<btb_infection_classes.size(); ++i)
    {
        btb_infection_classes[i] += n_in_btb_classes[i];
        current_size += n_in_btb_classes[i]; //Also keep track of the new total size.
    }

    //Select animals to remove from the classes to get the prem back to its original size
    //by making a sample from a multiv. hypergeometric distr.
    std::vector<unsigned int> removals;
    draw_multivariate_hypergeometric(current_size - prem_size, btb_infection_classes, removals);

    for(size_t i=0; i<btb_infection_classes.size(); ++i)
    {
        btb_infection_classes[i] -= removals[i];
    }
    btb_n_infectious = btb_infection_classes[4] + btb_infection_classes[5];
    verify_premsize();
}

void Prem_status::btb_subtract_from_inf_classes(size_t n, std::vector<unsigned int>& n_in_btb_classes)
{
    //Get the number
    size_t prem_size = this->get_spCountsCurrentQuarter().at(farm_type->get_species()); //"True" size.

    verify_premsize();
    n = std::min(n, prem_size); //Cant send more than whats currently on prem.
    draw_multivariate_hypergeometric(n, btb_infection_classes, n_in_btb_classes);
    //Update prevalence in this prem according to changes in animal numbers.
    //This is done by adding back as many animals as was removed to the classes
    //with probabilities of going into each separate class based on the new
    //prevalence in the classes. If there are no animals left, they are all
    //added to the sus class.

    int new_size = prem_size - n; //This many remains.
    if(new_size > 0)
    {
        //Update numbers by removing the shipped animals.
        for(size_t i=0; i<btb_infection_classes.size(); ++i)
        {
            btb_infection_classes[i] -= n_in_btb_classes[i];
        }
        //And replace with susc. animals
        btb_infection_classes[0] += n;
    }
    else
    {
         //All animals were sent off, replace all with sus animals.
         btb_infection_classes.clear();
         btb_infection_classes.resize(6, 0);
         btb_infection_classes[0] = n;
    }
    btb_n_infectious = btb_infection_classes[4] + btb_infection_classes[5];
    verify_premsize();
}

int Prem_status::btb_subtract_slaughter_shipment(size_t n_to_slaughter)
{
    verify_premsize();
    std::vector<unsigned int> to_slaughter(6, 0);
    btb_subtract_from_inf_classes(n_to_slaughter, to_slaughter);
    verify_premsize();
    return std::accumulate(to_slaughter.begin()+1, to_slaughter.end(), 0);
}

void Prem_status::verify_premsize()
{
    size_t prem_size = this->get_spCountsCurrentQuarter().at(farm_type->get_species()); //"True" size.
    size_t current_size = std::accumulate(btb_infection_classes.begin(), btb_infection_classes.end(), 0);
    if(prem_size != current_size)
    {
        std::cout << "Sum of individual classes is different from prem size. "
                  << prem_size << ", " << current_size << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Prem_status::btb_remove_all_infected()
{
    verify_premsize();
    for(size_t i=1; i<btb_infection_classes.size(); ++i)
    {
        btb_infection_classes[i] = 0;
    }
    btb_infection_classes[0] = this->get_spCountsCurrentQuarter().at(farm_type->get_species());
    btb_n_infectious = 0;
    verify_premsize();
}

bool Prem_status::btb_update_WH_spread(int quarter_idx)
{
    double tstep_frac = 1.0 / 12.0;

    //Some aliases for readability.
    unsigned int& nS = btb_infection_classes[0];
    unsigned int& nE1U = btb_infection_classes[1];
    unsigned int& nE2U = btb_infection_classes[2];
    unsigned int& nE2R = btb_infection_classes[3];
    unsigned int& nIU = btb_infection_classes[4];
    unsigned int& nIR = btb_infection_classes[5];

    unsigned int& nI = btb_n_infectious;
    nI = nIU + nIR;
    unsigned int n = nS + nE1U + nE2U + nE2R + nI;
    btb_n_new_WH_infections = 0;
    //Since the WH function is called at the beginning of every time step, here is one place
    //where the following variable can be reset before external infections happens. It's unrelated
    //to the WH spread process, though. Also adds the number of seeded animals since last timestep.
    btb_n_new_external_infections = 0;

    if(nS == n)
    {
        return false; //No need for WH-spread if only susceptibles.
    }
    else if(this->is_market() and parameters->within_spread_at_markets == false)
    {
        //One assumption we've made is that there is no bTB transmission going
        //on within markets. For bTB markets act only as redistributors of
        //shipments. But should any exposed animals that enter a market still
        //go through the entire transmission process as usual? Deaths? Births?
        //Or should any exposed animals that enter a market become 'frozen'
        //with regard to the development through the infection stages until
        //they are shipped off to the next premises?
        return true; //Returning true here means that infected animals at markets
        //are on hold until sent off.
    }

    //////////////////////////
    double d_n = double(n);
    double premsize = double(this->get_spCountsCurrentQuarter().at(farm_type->get_species()));
    verify_premsize();

    //////////////////////////////////
    // Pop. dynamics //
    //////////////////////////////////

    //My reasoning:
    //Import and export has been removed and replaced by USAMM shipments and shipments to slaughter.
    //In order to keep constant size, the number of newborn animals is assumed equal to the number
    //of dead animals. The animals that die due to mortality and replacement by newborn are subtracted
    //from all classes and replaced by susceptibles.
    //Animals from all classes are assumed to contribute to births, and they contribute before being
    //removed (so an animal is assumed to give birth, then die). Newborn go into susceptible,
    //but in order to keep the size constant, the same number as newborn need to be removed
    //from all classes. These replacements are made on the population numbers *after* the exported/dead
    //animals have been subtracted, not sure if this is the correct way to go.



    //For the transition periods (gamma distributed) we could perhaps use the mean and get rid
    //of the individual variation? It's still a stochastic process since the transition
    //period informs a rate of transitioning.

    //Additions from births. Only susceptibles.
    double birth_rate = 0.0; //Omega
    double removal_rate = 4.0 * btb_p->export_rate_v[quarter_idx] * (d_n / premsize) + btb_p->beef_mortality_rate;
    if(farm_type->get_species() == "dairy")
    {
        birth_rate = std::max((btb_p->dairy_birth_rate + btb_p->export_rate + btb_p->dairy_mortality_rate - btb_p->dairy_import_rate) *
                              (premsize / d_n), 0.0);
        removal_rate = 4.0 * btb_p->export_rate_v[quarter_idx] * (d_n / premsize) + btb_p->dairy_mortality_rate;
    }
    else if(quarter_idx == 1 and farm_type->get_species() == "beef") // Apr-Jun
    {
        birth_rate = std::max((btb_p->export_rate + btb_p->beef_mortality_rate - btb_p->beef_import_rate) *
                              (premsize / d_n), 0.0) * 4.0;
    }

    double birth_prob = 1.0 - std::exp(-birth_rate*tstep_frac);
    int n_newborn = draw_binom(n, birth_prob);

    //Removals due to export and death equally from all classes.
    double removal_p = 1.0 - std::exp(-removal_rate*tstep_frac);
    //Number of susc removed, not really of interest, since they are replaced by new susc anyway.
    unsigned int rem_E1U = draw_binom(nE1U, removal_p);
    unsigned int rem_E2U = draw_binom(nE2U, removal_p);
    unsigned int rem_E2R = draw_binom(nE2R, removal_p);
    unsigned int rem_IU = draw_binom(nIU, removal_p);
    unsigned int rem_IR = draw_binom(nIR, removal_p);



    //Now remove all exported/dead animals and add replacements to susceptibles.
    nS += rem_E1U + rem_E2U + rem_E2R + rem_IU + rem_IR;
    nE1U += -rem_E1U;
    nE2U += -rem_E2U;
    nE2R += -rem_E2R;
    nIU += -rem_IU;
    nIR += -rem_IR;
    verify_premsize();
    //With the updated sizes, replace animals with newborn with equal animal level probability
    //(the probability of replacing an animal of a specific class is proportional to the proportion
    //of the total population that is in that class).
    std::vector<unsigned int> to_remove;
    draw_multivariate_hypergeometric(n_newborn, {nS, nE1U, nE2U, nE2R, nIU, nIR}, to_remove);
    //We don't need to care about element 0 which is the number of susc replaced by newborn, since
    //that's equivalent of susc being replaced by susc.
    rem_E1U = to_remove[1];
    rem_E2U = to_remove[2];
    rem_E2R = to_remove[3];
    rem_IU = to_remove[4];
    rem_IR = to_remove[5];

    nS += rem_E1U + rem_E2U + rem_E2R + rem_IU + rem_IR;
    nE1U += -rem_E1U;
    nE2U += -rem_E2U;
    nE2R += -rem_E2R;
    nIU += -rem_IU;
    nIR += -rem_IR;

    nI = nIU + nIR;
    n = nS + nE1U + nE2U + nE2R + nI;
    d_n = double(n);
    verify_premsize();

    ////////////////////////////
    // Infections //
    ////////////////////////////

    //Calculate the individual probability for each susceptible animal to become infected.
    double infection_prob = 0.0;
    if(n > 0)
    {
        double inf_rate = btb_p->transmission_rate * tstep_frac *
                          (btb_p->cattle_cattle_contact_rate * nI / d_n +
                           btb_p->contact_rate_v[quarter_idx] +
                           btb_p->cattle_fomite_contact_rate);
        infection_prob = 1.0 - std::exp(-inf_rate);
    }

    //How many gets infected.
    btb_n_new_WH_infections = draw_binom(nS, infection_prob);

    //Animals transitioning E1.
    unsigned int n_trans_from_E1U = 0;
    unsigned int n_trans_from_E1U_to_E2R = 0;
    unsigned int n_trans_from_E1U_to_E2U = 0;
    if(nE1U > 0)
    {
        for(size_t i=0; i<nE1U; ++i)
        {
            double exp_period = draw_gamma(btb_p->sigma_1_shape, btb_p->sigma_1_scale);
            double leave_prob = 1.0 - std::exp(-exp_period*tstep_frac);
            if(uniform_rand() < leave_prob)
            {
                ++n_trans_from_E1U;
            }
        }
    }

    //Split them between E2U and E2R.
    if(n_trans_from_E1U > 0)
    {
        n_trans_from_E1U_to_E2R = draw_binom(n_trans_from_E1U,
                                             btb_p->test_sensitivity * btb_p->response_prob_1);
        n_trans_from_E1U_to_E2U = n_trans_from_E1U - n_trans_from_E1U_to_E2R;
    }

    //Animals transitioning from E2U
    unsigned int n_trans_from_E2U = 0;
    unsigned int n_trans_from_E2U_to_IU = 0;
    unsigned int n_trans_from_E2U_to_E2R = 0;
    if(nE2U > 0)
    {
        std::vector<double> exp_periods(nE2U, 0.0);
        std::vector<double> unr_periods(nE2U, 0.0);
        std::vector<int> which_transitions(nE2U, 0);

        for(size_t i=0; i<nE2U; ++i)
        {
            //Simulate exposed period for each animal.
            double exp_period = draw_gamma(btb_p->sigma_2_star_shape, btb_p->sigma_2_star_scale);
            exp_periods[i] = exp_period;
            //Simulate a period of being unreactive for each animal.
            double unr_period = draw_gamma(btb_p->delta_1_shape, btb_p->delta_1_scale);
            unr_periods[i] = unr_period;
            //Probability to become reactive or infectious, i.e. leaving E2U (according to Amanda, transitioning to infectious AND reactive is not allowed)
            double p_for_R_or_I = 1.0 - std::exp(-(exp_period * tstep_frac +
                                                   unr_period * tstep_frac * btb_p->test_sensitivity *
                                                   btb_p->response_prob_2));
            if(uniform_rand() < p_for_R_or_I)
            {
                ++n_trans_from_E2U;
                which_transitions[i] = 1;
            }
        }

        if(n_trans_from_E2U > 0)
        {
            for(size_t i=0; i<nE2U; ++i)
            {
                if(which_transitions[i] == 1)
                {
                    double p_for_I = exp_periods[i] /
                                     (exp_periods[i] + unr_periods[i] * btb_p->test_sensitivity * btb_p->response_prob_2);
                    if(uniform_rand() < p_for_I)
                    {
                        ++n_trans_from_E2U_to_IU;
                    }
                }
            }
            n_trans_from_E2U_to_E2R = n_trans_from_E2U - n_trans_from_E2U_to_IU;
        }
    }

    //Animals transitioning from E2R to IR
    unsigned int n_trans_from_E2R_to_IR = 0;
    if(nE2R > 0)
    {
        for(size_t i=0; i<nE2R; ++i)
        {
            //Simulate exposed period for each animal
            double exp_period = draw_gamma(btb_p->sigma_2_shape, btb_p->sigma_2_scale);
            double leave_E2R_prob = 1.0 - std::exp(-exp_period*tstep_frac);
            if(uniform_rand() < leave_E2R_prob)
            {
                ++n_trans_from_E2R_to_IR;
            }
        }
    }

    //Animals transitioning from IU to IR
    unsigned int n_trans_from_IU_to_IR = 0;
    //Simulate exposed periods.
    if(nIU > 0)
    {
        for(size_t i=0; i<nIU; ++i)
        {
            double exp_period = draw_gamma(btb_p->delta_2_shape, btb_p->delta_2_scale);
            double leave_IU_prob = 1.0 - std::exp(-exp_period * btb_p->test_sensitivity *
                                                  btb_p->response_prob_3 * tstep_frac);
            if(uniform_rand() < leave_IU_prob)
            {
                ++n_trans_from_IU_to_IR;
            }
        }
    }

    //Update all classes' sizes
    nS +=  -btb_n_new_WH_infections;
    nE1U += btb_n_new_WH_infections - n_trans_from_E1U;
    nE2U += n_trans_from_E1U_to_E2U - n_trans_from_E2U;
    nE2R += n_trans_from_E1U_to_E2R + n_trans_from_E2U_to_E2R - n_trans_from_E2R_to_IR;
    nIU += n_trans_from_E2U_to_IU - n_trans_from_IU_to_IR;
    nIR += n_trans_from_E2R_to_IR + n_trans_from_IU_to_IR;
    nI = nIU + nIR;
    n = nS + nE1U + nE2U + nE2R + nIU + nIR;

    //New true negs. and false positives. This should be moved to the code where the test is performed. If it remains here, it's not updated after changes to the numbers in the classes, such as when receiving a shipment for instance.
//    btb_true_neg = draw_binom(nS, btb_p->test_specificity * btb_p->specificity_scaling_factor);
//    btb_false_pos = nS - btb_true_neg;

    verify_premsize();
    if(nS == n)
    {
        return false; //No longer infected.
    }
    return true; //Still infected.
}



// This is a backup of the c++ implementation of the WH model for btb,
// before I removed everything that has anything to do with pop dynamics
// that are handled by quarterly flaps estimates and shipments instead.
// /S
//void Prem_status::btb_update_WH_spread(int quarter_idx)
//{
//    double tstep_frac = 1.0 / 12.0;
//
//    //Some aliases for readability.
//    unsigned int& nS = btb_infection_classes[0];
//    unsigned int& nE1U = btb_infection_classes[1];
//    unsigned int& nE2U = btb_infection_classes[2];
//    unsigned int& nE2R = btb_infection_classes[3];
//    unsigned int& nIU = btb_infection_classes[4];
//    unsigned int& nIR = btb_infection_classes[5];
//
//    unsigned int& nI = btb_n_infected;
//    nI = nIU + nIR;
//    unsigned int n = nS + nE1U + nE2U + nE2R + nI;
//
//    //////////////////////////
//    double d_n = double(n);
//    double premsize = double(speciesCounts[farm_type->get_species()]);
//
//
//    //////////////////////////////////
//    // Pop. dynamics //
//    //////////////////////////////////
//
//    double birth_rate = 0.0; //Omega
//    double removal_rate = 4.0 * btb_p->export_rate_v[quarter_idx] * (d_n / premsize) + btb_p->mortality_rate;
//    if(farm_type->get_species() == "dairy")
//    {
//        birth_rate = std::max((btb_p->dairy_birth_rate + btb_p->export_rate + btb_p->mortality_rate - btb_p->import_rate) *
//                              (premsize / d_n), 0.0);
//    }
//    else if(quarter_idx == 1 and farm_type->get_species() == "beef") // Apr-Jun
//    {
//        birth_rate = std::max((btb_p->export_rate + btb_p->mortality_rate - btb_p->import_rate) *
//                              (premsize / d_n), 0.0) * 4.0;
//    }
//
//    //Removals due to export and death equally from all classes.
//    double removal_p = 1.0 - std::exp(-removal_rate*tstep_frac);
//    unsigned int rem_S = draw_binom(nS, removal_p);
//    unsigned int rem_E1U = draw_binom(nE1U, removal_p);
//    unsigned int rem_E2U = draw_binom(nE2U, removal_p);
//    unsigned int rem_E2R = draw_binom(nE2R, removal_p);
//    unsigned int rem_IU = draw_binom(nIU, removal_p);
//    unsigned int rem_IR = draw_binom(nIR, removal_p);
//
//    //Additions from import.
//    double import_prob = 1.0 - std::exp(-btb_p->import_rate * tstep_frac);
//    unsigned int n_cur = ( (n == 0) ? premsize : n );
//    unsigned int imported = draw_binom(n_cur, import_prob);
//    //Number of imported animals that are susceptible.
//    unsigned int imp_S = draw_binom(imported, 1.0 - btb_p->prop_infected_imports);
//    //Number of imported infected animals that are reactive
//    unsigned int imp_reactive = draw_binom(imported - imp_S, btb_p->test_sensitivity);
//    //Out of which some are infectious...
//    unsigned int imp_IR = draw_binom(imp_reactive, btb_p->prop_reactive_are_infectious);
//    //...and the remaining are exposed.
//    unsigned int imp_E2R = imp_reactive - imp_IR;
//    //Number of imported to unreactive.
//    unsigned int imp_unreactive = imported - (imp_S + imp_reactive);
//    //Of which some are infectious...
//    unsigned int imp_IU = draw_binom(imp_unreactive, btb_p->prop_unreactive_are_infectious);
//    //...and exposed.
//    unsigned int imp_EU = imp_unreactive - imp_IU;
//    //Of which some have an immune response
//    unsigned int imp_E1U = draw_binom(imp_EU, btb_p->prop_unreactive_exposed_wo_immune_resp);
//    //...and some don't.
//    unsigned int imp_E2U = imp_EU - imp_E1U;
//
//    //Additions from births. Only susceptibles.
//    double birth_prob = 1.0 - std::exp(-birth_rate*tstep_frac);
//    unsigned int births_S = draw_binom(n, birth_prob);
//
//    //Add up additions to classes
//    unsigned int add_S = imp_S + births_S;
//    unsigned int add_E1U = imp_E1U;
//    unsigned int add_E2U = imp_E2U;
//    unsigned int add_E2R = imp_E2R;
//    unsigned int add_IU = imp_IU;
//    unsigned int add_IR = imp_IR;
//
//
//    ////////////////////////////
//    // Infections //
//    ////////////////////////////
//
//    //Calculate the individual probability for each susceptible animal to become infected.
//    double infection_prob = 0.0;
//    if(n > 0)
//    {
//        double inf_rate = btb_p->transmission_rate * tstep_frac * (btb_p->cattle_cattle_contact_rate * nI / d_n +
//                                                       btb_p->contact_rate_v[quarter_idx] +
//                                                       btb_p->cattle_fomite_contact_rate);
//        infection_prob = 1.0 - std::exp(-inf_rate);
//    }
//
//    //How many gets infected. *Is this wrong? Shouldn't it be binom(1,nS-m_S, b_no), otherwise you'd get more infectious than expected left after removing animals.
//    unsigned int n_newly_infected = 0;
//    n_newly_infected = draw_binom(nS, infection_prob);
//    n_newly_infected = std::max(0, int(std::min(nS - rem_S, n_newly_infected)));
//
//    //Animals transitioning E1.
//    unsigned int n_trans_from_E1U = 0;
//    unsigned int n_trans_from_E1U_to_E2R = 0;
//    unsigned int n_trans_from_E1U_to_E2U = 0;
//    if(nE1U > 0)
//    {
//    //            std::vector<double> exp_periods;
//    //            s1_gamma.draw_n(nE1U, exp_periods);
//
//        for(int i=0; i<nE1U; ++i)
//        {
//    //                double exp_period = draw_gamma(btb_p->sigma_1_shape, btb_p->sigma_1_scale);
//            double exp_period = draw_gamma(btb_p->sigma_1_shape, btb_p->sigma_1_scale);
//    //                double exp_period = exp_periods[i];
//            double leave_prob = 1.0 - std::exp(-exp_period*tstep_frac);
//            if(uniform_rand() < leave_prob)
//            {
//                ++n_trans_from_E1U;
//            }
//        }
//    }
//    n_trans_from_E1U = std::max(0, int(std::min(nE1U - rem_E1U, n_trans_from_E1U))); //Wrong again?
//
//    //Split them between E2U and E2R.
//    if(n_trans_from_E1U > 0)
//    {
//        n_trans_from_E1U_to_E2R = draw_binom(n_trans_from_E1U,
//                                             btb_p->test_sensitivity * btb_p->response_prob_1);
//        n_trans_from_E1U_to_E2U = n_trans_from_E1U - n_trans_from_E1U_to_E2R;
//    }
//
//    //Animals transitioning from E2U
//    unsigned int n_trans_from_E2U = 0;
//    unsigned int n_trans_from_E2U_to_IU = 0;
//    unsigned int n_trans_from_E2U_to_E2R = 0;
//    if(nE2U > 0)
//    {
//        std::vector<double> exp_periods(nE2U, 0.0);
//        std::vector<double> unr_periods(nE2U, 0.0);
//        std::vector<int> which_transitions(nE2U, 0);
//
//    //            s2star_gamma.draw_n(nE2U, exp_periods);
//    //            d1_gamma.draw_n(nE2U, unr_periods);
//        for(int i=0; i<nE2U; ++i)
//        {
//            //Simulate exposed period for each animal.
//    //                double exp_period = draw_gamma(btb_p->sigma_2_star_shape, btb_p->sigma_2_star_scale);
//            double exp_period = draw_gamma(btb_p->sigma_2_star_shape, btb_p->sigma_2_star_scale);
//            exp_periods[i] = exp_period;
//    //                double exp_period = exp_periods[i];
//            //Simulate a period of being unreactive for each animal.
//    //                double unr_period = draw_gamma(btb_p->delta_1_shape, btb_p->delta_1_scale);
//            double unr_period = draw_gamma(btb_p->delta_1_shape, btb_p->delta_1_scale);
//            unr_periods[i] = unr_period;
//    //                double unr_period = unr_periods[i];
//            //Probability to become reactive or infectious, i.e. leaving E2U (according to Amanda, transitioning to infectious AND reactive is not allowed)
//            double p_for_R_or_I = 1.0 - std::exp(-(exp_period * tstep_frac +
//                                                   unr_period * tstep_frac * btb_p->test_sensitivity *
//                                                   btb_p->response_prob_2));
//            if(uniform_rand() < p_for_R_or_I)
//            {
//                ++n_trans_from_E2U;
//                which_transitions[i] = 1;
//            }
//        }
//
//        n_trans_from_E2U = std::max(0, int(std::min(nE2U - rem_E2U, n_trans_from_E2U))); //Wrong again...
//
//        if(n_trans_from_E2U > 0)
//        {
//            for(int i=0; i<nE2U; ++i)
//            {
//                if(which_transitions[i] == 1)
//                {
//                    double p_for_I = exp_periods[i] /
//                                     (exp_periods[i] + unr_periods[i] * btb_p->test_sensitivity * btb_p->response_prob_2);
//                    if(uniform_rand() < p_for_I)
//                    {
//                        ++n_trans_from_E2U_to_IU;
//                    }
//                }
//            }
//            n_trans_from_E2U_to_E2R = n_trans_from_E2U - n_trans_from_E2U_to_IU;
//        }
//    }
//
//    //Animals transitioning from E2R to IR
//    unsigned int n_trans_from_E2R_to_IR = 0;
//    if(nE2R > 0)
//    {
//    //            std::vector<double> exp_periods;
//    //            s2_gamma.draw_n(nE2R, exp_periods);
//        for(int i=0; i<nE2R; ++i)
//        {
//            //Simulate exposed period for each animal
//    //                double exp_period = draw_gamma(btb_p->sigma_2_shape, btb_p->sigma_2_scale);
//            double exp_period = draw_gamma(btb_p->sigma_2_shape, btb_p->sigma_2_scale);
//    //                double exp_period = exp_periods[i];
//            double leave_E2R_prob = 1.0 - std::exp(-exp_period*tstep_frac);
//            if(uniform_rand() < leave_E2R_prob)
//            {
//                ++n_trans_from_E2R_to_IR;
//            }
//        }
//    }
//    n_trans_from_E2R_to_IR = std::max(0, int(std::min(nE2R - rem_E2R, n_trans_from_E2R_to_IR))); //Wrong...
//
//    //Animals transitioning from IU to IR
//    unsigned int n_trans_from_IU_to_IR = 0;
//    //Simulate exposed periods.
//    if(nIU > 0)
//    {
//    //            std::vector<double> exp_periods;
//    //            d2_gamma.draw_n(nIU, exp_periods);
//        for(int i=0; i<nIU; ++i)
//        {
//    //                double exp_period = draw_gamma(btb_p->delta_2_shape, btb_p->delta_2_scale);
//            double exp_period = draw_gamma(btb_p->delta_2_shape, btb_p->delta_2_scale);
//    //                double exp_period = exp_periods[i];
//            double leave_IU_prob = 1.0 - std::exp(-exp_period * btb_p->test_sensitivity *
//                                                  btb_p->response_prob_3 * tstep_frac);
//            if(uniform_rand() < leave_IU_prob)
//            {
//                ++n_trans_from_IU_to_IR;
//            }
//        }
//    }
//    n_trans_from_IU_to_IR  = std::max(0, int(std::min(nIU - rem_IU, n_trans_from_IU_to_IR ))); //Wrong...
//
//    //Update all classes' sizes
//    nS += add_S - n_newly_infected - rem_S;
//    nE1U += add_E1U + n_newly_infected - n_trans_from_E1U - rem_E1U;
//    nE2U += add_E2U + n_trans_from_E1U_to_E2U - n_trans_from_E2U - rem_E2U;
//    nE2R += add_E2R + n_trans_from_E1U_to_E2R + n_trans_from_E2U_to_E2R - n_trans_from_E2R_to_IR - rem_E2R;
//    nIU += add_IU + n_trans_from_E2U_to_IU - n_trans_from_IU_to_IR - rem_IU;
//    nIR += add_IR + n_trans_from_E2R_to_IR + n_trans_from_IU_to_IR - rem_IR;
//    nI = nIU + nIR;
//    n = nS + nE1U + nE2U + nE2R + nIU + nIR;
//
//    //New true negs. and false positives
//    btb_true_neg = draw_binom(nS, btb_p->test_specificity * btb_p->specificity_scaling_factor);
//    btb_false_pos = nS - btb_true_neg;
//}
