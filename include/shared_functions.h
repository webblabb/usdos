//© 2024 Colorado State University
/// \file

#ifndef shared_functions_h
#define shared_functions_h

#include <algorithm>
#include <random> // for random number generator
#include <chrono> // for random number generator
#include <cmath> // for std::sqrt in gKernel, floor in randomFrom
#include <fstream> // for printing
#include <iostream> // for troubleshooting output
#include <sstream>
#include <unordered_map>
#include <gsl/gsl_rng.h>
#include "gsl/gsl_randist.h"


#include "Farm.h"

    //Typedefs of some multidimensional vectors.
    typedef std::vector<double> Vec_d_1d;
    typedef std::vector<Vec_d_1d> Vec_d_2d;
    typedef std::vector<Vec_d_2d> Vec_d_3d;
    typedef std::vector<Vec_d_3d> Vec_d_4d;
    typedef std::vector<Vec_d_4d> Vec_d_5d;

    typedef std::vector<int> Vec_i_1d;
    typedef std::vector<Vec_i_1d> Vec_i_2d;
    typedef std::vector<Vec_i_2d> Vec_i_3d;

    typedef std::map<std::string, std::vector<double>> str_vec_map;

    //Shared functions
	double uniform_rand(); ///< Uniform distribution random number generator
	int rand_int(int lo, int hi); ///< Uniform integer distribution rng.
	double draw_norm(double mean, double stdev); ///< Draw a random number from a normal distribution with specified mean and standard deviation.
	int draw_binom(int, double); ///< Draw number of successes from a binomial distribution
	int draw_poisson(double lambda); ///<Generate a random number from a poisson dist. with given rate.
	double draw_gamma(double shape, double scale); ///<Generate a random number from a gamma distribution with parameters shape and scale.
	void draw_multinomial(int n, const std::vector<double>& weights, std::vector<unsigned int>& outcome); ///<Generate a random sample from a multinomial distribution and save to outcome.
	void draw_multivariate_hypergeometric(int n, const std::vector<unsigned int>& bins, std::vector<unsigned int>& outcome); ///<Generate a sample from a multivariate hypergeometric distribution with a vector describing the number of objects in each bin. Same as multinomial but without replacement.
	double draw_beta(double alpha, double beta); ///<Generate a random number from a beta distribution.
	unsigned int generate_distribution_seed(); ///<Generates a number that can be used in a random number generator based on the current time.
	size_t get_day_of_year(size_t current_timestep, size_t start_day); ///<Given the current time step and the start day of the simulation, returns the current day of the year.
	size_t get_month_of_year(size_t day_of_year);
	int get_quarter(size_t current_day_of_year);
	double oneMinusExp(double); ///< Calculates \f$1 - e^x\f$ using a two-term Taylor approximation for x<1e-5
 	int normDelay(std::tuple<double, double>&); ///< Return a period of time, drawn from a normal distribution
	std::vector<std::string>
		split(const std::string&, char, std::vector<std::string>&);
	std::vector<std::string>
		split(const std::string&, char);
    void skipBOM(std::ifstream &in);
	std::string to_string(Farm*);
	std::vector<std::string> semicolonStringToStringVec(std::string&); ///< Converts semicolon-separated string to vector of strings
	std::vector<double> stringToNumVec(std::string&); ///< Converts comma-separated string to vector of doubles
	std::vector<int> stringToIntVec(std::string&); ///< Converts comma-separated string to vector of integers
	std::vector<std::string> stringToStringVec(std::string&); ///< Converts comma-separated string to vector of strings
	std::string vecToCommaSepString(const std::vector<int>); ///< Converts vector of integers to a comma-separated string
	std::string vecToCommaSepString(const std::vector<std::string>); ///< Overloaded version converts vector of strings to a comma-separated string
	void addItemTab(std::string&, int); ///< Adds tab after an integer (converted to character)
	void addItemTab(std::string&, double); ///< Overloaded version adds tab after a double (converted to character)
	void addItemTab(std::string&, std::string); ///< Overloaded version adds tab after a string
	void printLine(std::string&, std::string&); ///< Generic print function used by a variety of output files
	unsigned int get_n_lines(std::ifstream& f); ///< Counts and returns the number of lines in a file.
    int read_table(std::string fname, char delim,bool header,
                   std::vector<std::vector<std::string>>& out,
                   char comment_char='#');
	std::string trim(const std::string &s); //Removes whitespace around string.
	size_t pairing_function(int x, int y); //Takes two natural numbers and returns a new number that is unique for that combination. Used for turning a combination of two prem classes into a single key for lookup in maps.


template<typename T>
T stringToNum(const std::string& text)
{
	std::istringstream ss(text);
	T result;
	if (! (ss>>result)){result = -1;}
	return result;
}

///> Chooses a single random element from a vector
///	\param[in]	vec		Vector of values from which to choose
template<typename T>
T randomFrom(std::vector<T>& vec)
{
	int maxSize = vec.size();
	double rUnif = uniform_rand();
	int rIndex = (int) floor(rUnif*maxSize);

	if (rUnif == 1){rIndex = vec.size()-1;} // Return last element (rather than out-of-range vec[maxSize])

	return vec[rIndex];
}

///> Chooses multiple random elements from a vector based on the Fisher-Yates shuffling algorithm.
/// num_random selected random values are copied to output,
/// selected values are swapped to the end of the vector so they're not selected again
/// Used in Grid_manager to select binomial-success farms
/// Used in Shipping_manager to select random premises in counties
/// \param[in]	elements	Vector of values from which to choose
///	\param[in]	num_random	Number of values to choose
/// \param[out] output1		Vector of randomly chosen values
template<typename T>
void random_unique(std::vector<T> elements, int num_random, std::vector<T>& output1)
	// elements not referenced (&) because we're rearranging it
{
	std::vector<T> output;
	output.reserve(num_random);
	int endIndex = elements.size();
	// endIndex separates non-selected values (elements [0, endIndex-1]) from selected values
	for (auto i = 1; i<= num_random; i++){
		// choose random number between 0 and 1
		double rUnif = uniform_rand();
		// scale up to endIndex, so r is an index in [0, endIndex)
		if (rUnif == 1 ){rUnif=0.999;} // avoids assigning actual endIndex value (out of range)
		int r = (int)floor(rUnif*endIndex);
		// copy value to output
		output.emplace_back(elements[r]);
		// swap r and endIndex-1 (last element of non-selected values)
		std::swap(elements[r], elements[endIndex-1]);
		endIndex--;
	}
 output.swap(output1);
}

///> Checks if an item is within a vector of items
template<typename T>
bool isWithin(const T target, const std::vector<T> vec)
{
	auto it = vec.begin();
	bool found = 0;
	while (it!=vec.end() && found==0){ // keep going while not done with vector and target not found
		if(*it == target){found = 1;}
		it++;
	}
	return found;
}

///> Determine which element's range a number falls into
/// Used in determining which (shipment) method to use at a given time
/// \param[in]	toMatch			Value to be matched
/// \param[in]	elementMaxes	Sorted vector of maximums for each element
/// \returns largest element of elementMaxes that is less than or = toMatch
template<typename T>
int whichElement(T& toMatch, std::vector<T>& elementMaxes)
{
	int match = -1; // the element that will be returned
	if (toMatch > elementMaxes.back()){
		std::cout<<"Error: (whichElement): value to match exceeds largest of comparison values. Exiting..."
		<< std::endl;
		exit(EXIT_FAILURE);}
	if (elementMaxes.size() < 1){
		std::cout << "Error (whichElement): Vector of element sizes < 1. Exiting..."
		<< std::endl;
		exit(EXIT_FAILURE);}
	if (elementMaxes.size()==1){match=0;}
	else{
		bool found = 0;
		unsigned int it = 1;
		while (it!=elementMaxes.size() && found == 0){
			if (toMatch>=elementMaxes[it-1] && toMatch<elementMaxes[it]){ // >= than previous, < current
				match = it-1; // subtract one to get the element below
				found = 1;
			}
			it++;
		}
		if (it==elementMaxes.size() && found == 0){
			std::cout << "Warning: (whichElement): Match not found.";}
	}
	return match;
}

///> Function that compares Farm IDs for sorting
/// Used with grid_cell*s in grid checker: stepThroughCells. Must be defined outside of
/// class in order to be used with std::sort
template<typename T>
inline bool sortByID(const T item1, const T item2)
{
	return (item1 -> get_id()) < (item2 -> get_id());
}

#endif
