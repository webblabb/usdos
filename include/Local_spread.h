#ifndef Local_spread_h
#define Local_spread_h

#include <iterator> // std::prev
#include <map>
#include "shared_functions.h" // for split; contains cmath, fstream, iostream

extern int verboseLevel;

/// Defines relationship between distance and transmission risk.
/// May be defined as a function or read in as an external file. Constructors perform 
/// one-time initial calculations (i.e. operations on parameters or squaring distances) 
/// and store them with the instance.
class Local_spread
{
	private:
		int verbose; ///< Can be set to override global setting for console output
		int kType; ///< Kernel type
		std::vector<double> kp; ///< Kernel parameters
		std::string datafile; ///< File containing distances and associated probabilities
		std::map<double,double> distProb; ///< Map with key of distance (m) squared, value is probability	
		
	public:
		///> Constructs a kernel from an equation (determined by variable kType)
		Local_spread(int kernelType, 
			std::vector<double> kparams = std::vector<double>());
		///> Constructs a kernel from external file
		Local_spread(int kernelType, 
			std::string fname);
		~Local_spread();
		///> Calculates or matches kernel value according to form defined at construction
		double atDistSq(double);
};
	
#endif // Local_spread_h