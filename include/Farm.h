#ifndef FARM_H
#define FARM_H

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <map> // for distances to neighbors
#include <string>
#include <unordered_map>
#include <utility> // for std::iter_swap in rem_probPreventExposure
#include <vector>
#include "Point.h"

class County;
class State;
class Farm_type;
struct Parameters;

extern int verboseLevel;

/// Describes a premises - one of these objects is created for each premises row read in
/// from the premises file.
class Farm
{
	protected: // allows access from derived class Prem_status
		int id,	///< Unique integer identifier read from premises file
			cellID; ///< Integer identifier of grid_cell assigned to this premises during grid creation
		double x_coordinate, ///< x-coordinate from projected longitude (same units as local spread kernel)
			y_coordinate, ///< y-coordinate from projected latitude (same units as local spread kernel)
			sus, ///< Calculated total susceptibility of this premises
			inf; ///< Calculated total infectiousness of this premises
        Point position;
		County* parent_county;
		State* parent_state;
		Farm_type* farm_type;
		std::string fips; ///< County identifier (FIPS code)
		double neighborRadiusCalculated; ///< If neighbors have been calculated and stored as distancesNeighbors, the radius searched (else 0)
		std::unordered_map< std::string, int > speciesCounts; ///< Numbers of animals of each type, keyed by types
		std::multimap<double, Farm*> distancesNeighbors; ///< Stored calculations of distances to neighboring premises
		double oweight; ///< The origin weight of this farm in relation to all other farms within the same state of the same type.
		double normalized_oweight; ///< The origin weight of this farm in relation to all other farms within the same state of the same type, normalized so the sum of all farms in state = 1.0.
		double dweight; ///< The destination weight of this farm in relation to all other farms within the same state of the same type.

	public:
		Farm(int, double, double, std::string);
		~Farm();
		int get_id() const; // inlined
		int get_cellID() const; // inlined
		Farm_type* get_farm_type() const; //inlined
		double get_x() const; // inlined
		double get_y() const; // inlined
		double get_sus() const; // inlined
		double get_inf() const; // inlined
        double get_sus_max() const; // inlined
        double get_inf_max() const; // inlined
        double get_normalized_oweight() const;
		double get_unnormalized_oweight() const;
		double get_unnormalized_dweight() const;
 		std::string get_fips() const; // inlined
 		const std::unordered_map< std::string, int >& get_spCounts(); // inlined
 		int get_size(const std::string species) const;
 		int get_size_allSpecies() const;
		County* get_parent_county() const; //Inlined
 		State* get_parent_state() const; //Inlined
 		const std::multimap<double, Farm*>* get_distancesNeighbors(); //inlined
 		double get_neighborRadiusCalculated() const; //inlined
		void set_cellID(const int cellID);
		void set_farm_type(Farm_type* in_type);
 		void set_speciesCount(const std::string, int);
 		void set_sus(const double);
 		void set_inf(const double);
 		void set_normalized_oweight(const double in_oweight);
 		void set_unnormalized_oweight(const double in_oweight);
 		void set_unnormalized_dweight(const double in_dweight);
		void set_parent_county(County* in_county);
		void set_distancesNeighbors(std::multimap<double, Farm*>&);
		void set_neighborRadiusCalculated(const double);

};

inline int Farm::get_id() const
{
	return id;
}
inline int Farm::get_cellID() const
{
	return cellID;
}
inline Farm_type* Farm::get_farm_type() const
{
  return farm_type;
}
inline double Farm::get_x() const
{
	return x_coordinate;
}
inline double Farm::get_y() const
{
	return y_coordinate;
}
inline double Farm::get_sus() const
{
	return sus;
}
inline double Farm::get_inf() const
{
	return inf;
}
inline double Farm::get_sus_max() const
{
    return sus;
}

inline double Farm::get_inf_max() const
{
    return inf;
}
inline std::string Farm::get_fips() const
{
	return fips;
}
const inline std::unordered_map< std::string, int >& Farm::get_spCounts()
{
	return speciesCounts;
}
inline County* Farm::get_parent_county() const
{
  return parent_county;
}
const inline std::multimap<double, Farm*>* Farm::get_distancesNeighbors()
{
	return &distancesNeighbors;
}
inline double Farm::get_neighborRadiusCalculated() const
{
	return neighborRadiusCalculated;
}


class Farm_type
{
public:
    Farm_type(std::string herd, std::vector<std::string> in_species, unsigned int index);
    ~Farm_type();
    std::string get_species() const; //inlined
    unsigned int get_index() const; //inlined
private:
    unsigned int index; // an index value relative to all other Farm_types
    std::string herd; // a binary code (ie '01') unique to this Farm_type, used to indicate species composition
    std::string species; //ie dairy, beef, name for this farm type

};

inline std::string Farm_type::get_species() const
{
    return species;
}
inline unsigned int Farm_type::get_index() const
{
    return index;
}

/// Derived class inheriting from Farm, containing additional info on replicate-specific statuses
/// and deleted at the end of each replicate. One object of this type is created for each
/// premises per simulation when it has any kind of status change (i.e. reporting, exposed,
/// prophylactic vaccination).
class Prem_status: public Farm
{
	private:
		std::string fileStatus; /// Current file status, one of: notDangerousContact, dangerousContact, reported
		std::string diseaseStatus; /// Current disease status, one of: exp, inf, imm
		std::vector<std::string> controlStatuses; ///< Vector of all control-related statuses ever effective for this premises
		std::unordered_map<std::string, std::string> statusWhenWaitlisted; /// Record of premises status at the time of addition to control waitlist - the reason why control was targeted. Key is control type, value is status.

		std::unordered_map<std::string, int> start; /// Start times for each status. Assumes unique names among file, disease, and control statuses due to mapping by status name.
		std::unordered_map<std::string, int> end; /// End times for each status. Assumes unique names among file, disease, and control statuses due to mapping by status name.

		std::vector<double> probPreventExposure; /// Probability of exposure being prevented, given effective control. One element for each control type currently effective.
		std::vector<double> probPreventTransmission; /// Probability of transmission being prevented, given effective control. One element for each control type currently effective.
		std::vector<Farm*> dangerousContacts; /// Farms that are dangerousContacts (set when premises is reported)
		std::unordered_map<Farm*, std::unordered_map<std::string, bool>> potentialDCs; /// Map keyed by pointer to potentialDC, then by disease status, then bool if potentialDC is that status when this premises is reported

		// These vectors are "aligned" - each element corresponds to the same element in the others. Guessed this was easier and more flexible to search than tuple elements.
		std::vector<int> expTime;
		std::vector<Farm*> expSource;
		std::vector<int> expRoute;
		std::vector<std::string> expBlocked;

        bool isVaccinated;
        std::unordered_map< std::string, int > currentSize; ///< Numbers of animals of each type, keyed by types
        std::unordered_map< std::string, int > currentSizeUnvaccinated; ///< Numbers of animals of each type that are unvaccinated, keyed by types

        double get_inf_partial(int t, const Parameters* p,
                               const std::unordered_map<std::string, double>& normInf_map,
                               const std::unordered_map<std::string, int>& sp_counts);

	public:
		Prem_status(Farm*);
		~Prem_status();

		void add_controlStatus(const std::string); //inlined
 		void add_exposureSource(int, Farm*, int, std::string); // confirmed record of source of infection, route, block

 		void set_fileStatus(const std::string); //inlined
 		void set_diseaseStatus(const std::string); //inlined
 		void set_start(const std::string, const int); //inlined - set start time for status
 		void set_end(const std::string, const int); //inlined - set end time for status
 		void set_statusWhenWaitlisted(const std::string, const std::string);//inlined

		double get_probPreventExposure() const;
 		void add_probPreventExposure(double); // inlined
 		void rem_probPreventExposure(double);

		double get_probPreventTransmission() const;
 		void add_probPreventTransmission(double); // inlined
		void rem_probPreventTransmission(double);

        bool get_isVaccinated() { return isVaccinated; }
        void vaccinate(double efficacy); //Determines the number of animals that are still susceptible after vaccination as  N_s~Bin(N_tot, 1-efficacy)
 		void unvaccinate(); //Removes effect of vaccination in the herd and sets N_s = N_tot.
 		const std::unordered_map<std::string, int>& get_currentSizeUnvaccinated();

 		std::string get_diseaseStatus() const; //inlined
 		std::string get_fileStatus() const; //inlined
 		int get_start(std::string) const;
 		int get_end(std::string) const;
 		std::vector<Farm*> get_dangerousContacts() const; //inlined
 		std::vector<std::string> get_controlStatuses() const; //inlined
 		bool is_exposureSource(Farm*);
 		void add_potentialDCInfo(Farm*, const std::unordered_map<std::string, bool>&);//inlined
 		std::unordered_map<Farm*, std::unordered_map<std::string, bool>>* get_potentialDCs();//inlined
 		void add_DC(Farm* f); //inlined
 		bool was_dcWhenWaitlisted(const std::string);
 		bool is_onWaitlist(const std::string);

 		bool beenExposed() const; //inlined - not used?

        double when_infected(); //returns earliest exposure that was not blocked

        void set_currentSize(const std::string species, int sp_count); //inlined
        int get_currentSize(const std::string species);

        double get_inf_partial_as_unvaccinated(int t, const Parameters* p,
                                               const std::unordered_map<std::string, double>& normInf_map);
        double get_inf_partial_as_vaccinated(int t, const Parameters* p,
                                             const std::unordered_map<std::string, double>& normInf_map);
};

inline std::string Prem_status::get_fileStatus() const
{
	return fileStatus;
}
inline std::string Prem_status::get_diseaseStatus() const
{
	return diseaseStatus;
}
inline std::vector<Farm*> Prem_status::get_dangerousContacts() const
{
	return dangerousContacts;
}
inline std::vector<std::string> Prem_status::get_controlStatuses() const
{
	return controlStatuses;
}
inline void Prem_status::set_fileStatus(const std::string status)
{
	fileStatus = status;
}
inline void Prem_status::set_diseaseStatus(const std::string stat)
{
	diseaseStatus = stat;
}
inline void Prem_status::add_controlStatus(const std::string status)
{
	controlStatuses.emplace_back(status);
}
inline void Prem_status::set_start(const std::string status, const int t)
{
	start[status] = t;
}
inline void Prem_status::set_end(const std::string status, const int t)
{
	end[status] = t;
}
inline void Prem_status::set_statusWhenWaitlisted(const std::string c_type, const std::string status)
{
	statusWhenWaitlisted[c_type] = status;
}

inline void Prem_status::add_probPreventExposure(double eff)
{
	probPreventExposure.emplace_back(eff);
}
inline void Prem_status::add_probPreventTransmission(double eff)
{
	probPreventTransmission.emplace_back(eff);
}
/// Checks whether or not a premises has been exposed to infection by checking for the presence of an "exp" start time
inline bool Prem_status::beenExposed() const
{
	return start.count("exp")==1;
}
inline std::unordered_map<Farm*, std::unordered_map<std::string, bool>>* Prem_status::get_potentialDCs()
{
	return &potentialDCs;
}
inline void Prem_status::add_DC(Farm* f)
{
	dangerousContacts.emplace_back(f);
}
inline bool Prem_status::is_onWaitlist(const std::string c_type)
{
	return statusWhenWaitlisted.count(c_type)>0;
}
inline void Prem_status::set_currentSize(const std::string species, int sp_count)
{
    currentSize[species] = sp_count;
}

#endif //FARM_H
