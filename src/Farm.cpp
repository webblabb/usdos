#include <stdio.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "Farm.h"
#include "County.h"
#include "shared_functions.h"
#include "State.h"
#include "File_manager.h"

Farm::Farm(int in_id, double in_x, double in_y, std::string in_fips)
	:
	id(in_id),
	cellID(-1),
	x_coordinate(in_x),
	y_coordinate(in_y),
	position(in_x, in_y),
	fips(in_fips),
	neighborRadiusCalculated(0)
{
}

Farm::~Farm()
{
}

double Farm::get_normalized_oweight() const
{
    return normalized_oweight;
}

double Farm::get_unnormalized_oweight() const
{
    return oweight;
}

double Farm::get_unnormalized_dweight() const
{
    return dweight;
}


int Farm::get_size(const std::string species) const
{
	int count = 0;
	if (speciesCounts.count(species) != 0){count = speciesCounts.at(species);}

	return count;
}
int Farm::get_size_allSpecies() const
{
	int count = 0;
	for (auto& i:speciesCounts){
		count += i.second;
	}
  return count;
}

void Farm::set_cellID(const int in_cellID)
{
	cellID = in_cellID;
}

void Farm::set_farm_type(Farm_type* in_type)
{
    farm_type = in_type;
}

void Farm::set_speciesCount(const std::string species, int sp_count)
{
	speciesCounts[species] = sp_count;
}

void Farm::set_sus(const double in_sus)
{
	sus = in_sus;
}

void Farm::set_inf(const double in_inf)
{
	inf = in_inf;
}

void Farm::set_normalized_oweight(const double in_normalized_oweight)
{
    normalized_oweight = in_normalized_oweight;
}

void Farm::set_unnormalized_oweight(const double in_oweight)
{
    oweight = in_oweight;
}

void Farm::set_unnormalized_dweight(const double in_dweight)
{
    dweight = in_dweight;
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
	isVaccinated(false),
	currentSizeUnvaccinated(speciesCounts)
{
	probPreventExposure.emplace_back(0.0);
	probPreventTransmission.emplace_back(0.0);
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
    for(auto& sp_count_pair : speciesCounts)
    {
        if(sp_count_pair.second > 0)
        {
            int n_unvaxed = draw_binom(sp_count_pair.second, 1.0 - efficacy);
            currentSizeUnvaccinated[sp_count_pair.first] = n_unvaxed;
        }
    }
    isVaccinated = true;
}

void Prem_status::unvaccinate()
{
    currentSizeUnvaccinated = speciesCounts;
    isVaccinated = false;
}

///Returns the number of unvaccinated animals, if the herd is vaccinated, it's
///the part of the herd for which the vaccine was inefficient; if the herd isn't
///vaccinated it's all animals.
const std::unordered_map<std::string, int>& Prem_status::get_currentSizeUnvaccinated()
{
    return currentSizeUnvaccinated;
}

/// Stores info for determining if potential DCs will be dangerousContacts
void Prem_status::add_potentialDCInfo(Farm* pdc, const std::unordered_map<std::string, bool>& dcInfo)
{
	if (potentialDCs.count(pdc)==0){
		potentialDCs[pdc] = dcInfo;
		if(verboseLevel>1){std::cout<<"F in Prem_status add_potentialDCInfo first if statement"<<std::endl;}

	} else { // in case DC relationship already established between this pair of farms,
	         // isDC status of 'true' should override any 'false' evaluations
		if(verboseLevel>1){std::cout<<"F in Prem_status add_potentialDCInfo first if statement else"<<std::endl;}

if(verboseLevel>1){std::cout<<"DC relationship already established, new values: ";}
		std::unordered_map<std::string, bool> existingPDC = potentialDCs.at(pdc);
		for (auto& s:existingPDC){ //s.first is status (i.e. sus, exp); s.second is bool
			if (dcInfo.at(s.first)==1){
				potentialDCs.at(pdc).at(s.first) = 1;
				if(verboseLevel>1){std::cout<<"F in Prem_status add_potentialDCInfo for loop"<<std::endl;}

			}
if(verboseLevel>1){std::cout<<s.first<<"="<<potentialDCs.at(pdc).at(s.first)<<", ";}
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

/// Returns true if premises was not reported at time of addition to waitlist
bool Prem_status::was_dcWhenWaitlisted(const std::string c_type){
	return (statusWhenWaitlisted.at(c_type) == "dangerousContact" || statusWhenWaitlisted.at(c_type) == "exposed");
}

/// Returns the time of infection for a farm
double  Prem_status::when_infected()
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



int Prem_status::get_currentSize(const std::string species)
{
    int count = 0;
    if (currentSize.count(species) != 0){count = currentSize.at(species);}

    return count;
}

///calculates current infectiousness for farm and current time t using
///the entire population size.
double Prem_status::get_inf_partial_as_unvaccinated(int t, const Parameters* p,
                                                    const std::unordered_map<std::string, double>& normInf_map)
{
    return get_inf_partial(t, p, normInf_map, this->get_spCounts());
}

///calculates current infectiousness for farm and current time t assuming
///that the farm has been vaccinated and using the part of population size for which
///vaccination was unsuccessful.
double Prem_status::get_inf_partial_as_vaccinated(int t, const Parameters* p,
                                                  const std::unordered_map<std::string, double>& normInf_map)
{
    return get_inf_partial(t, p, normInf_map, this->get_currentSizeUnvaccinated());
}

///calculates current infectiousness for farm and current time t
double Prem_status::get_inf_partial(int t, const Parameters* p,
                                    const std::unordered_map<std::string, double>& normInf_map,
                                    const std::unordered_map<std::string, int>& sp_counts)
{
    //find time of infection
    double t0=this->when_infected();

    ////fixed parameters
    const std::vector<double>& partialParams = p->partialParams;
    double r0=partialParams.at(0);
    double r1= partialParams.at(1);
    double gamma=partialParams.at(2);;
    double tS0 =partialParams.at(3);
    double a=partialParams.at(4);
    double b=partialParams.at(5);

    double I1=0.0;
    double I2=0.0;

    //tsigma is the same as the latency duration
    double tsigma=std::get<0>(p->latencyParams);

    //duration of time infected
    //must be converted to double for below calculation
    double tInf= (double(t)-t0);

    double premInf = 0.0;

    // for each species type, calcuate the number of infectious animals
    // loop through species on prems
    for (auto& sp_N_pair : sp_counts){

        double NInf=0.0; //comment back in

        double N = double(sp_N_pair.second); // i.e. get_size("beef") gets # of beef cattle on premises

        ////calculate growith rate
        double r= N*(r0*(tInf-1)+r1*std::pow(tInf,2));

        ////evaluate growth curve

        if(tInf<tsigma) {NInf=I1;}

        if((tInf>=tsigma)&&tInf<(tsigma+tS0)){ NInf += (I1+I2*(tInf-tsigma)/tS0+(r+gamma*r*(tS0-(tInf-tsigma)))/(tS0*std::pow(gamma,2))-((r+gamma*r*tS0)/(tS0*std::pow(gamma,2)))*exp(-gamma*(tInf-tsigma)));}

        if(tInf>=(tsigma+tS0)){ NInf += (I1+I2+(r*std::exp(gamma*tS0)-r-gamma*r*tS0)/(tS0*std::pow(gamma,2))*std::exp(-gamma*(tInf-tsigma)));}

        double normInf_val=normInf_map.at(sp_N_pair.first);

        //scale up transmissibility norminf
        double factor=a+b*N;

        double spInf =(factor*normInf_val)*std::pow(NInf, p->infExponents.at(sp_N_pair.first)); // susceptibility value for this species/type

        premInf += spInf; // add this species to the total for this premises

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
    return premInf;
}
