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
		std::vector<double> kp_q1; ///< Kernel parameters for quarter 1
		std::vector<double> kp_q2; ///< Kernel parameters for quarter 2
		std::vector<double> kp_q3; ///< Kernel parameters for quarter 3
		std::vector<double> kp_q4; ///< Kernel parameters for quarter 4
		std::string datafile; ///< File containing distances and associated probabilities
		std::map<double,double> distProb; ///< Map with key of distance (m) squared, value is probability

		std::vector<double> V_btb_local, V_btb_wildlife; ///< Stores volumes of original btb kernels so they can be normalized to have vol = 1.0, which is necessary for a mixture kernel with local and wildlife trasmission.
		double btb_mixed_kernel(double d, double wildl_dens, double lp1, double lp2,
                                double wp1, double wp2, double x, int quarter_idx); //Mix of local (L) and wildlife (W), x controls proportion of entire kernel that is wildlife. M(d) = x*W(d) + (1-x)*L(d), automatically normalized so that volume of M == volume of L on it's own.
		static double btb_local(double d, double phi, double alpha);
		static double btb_local_shell(double d, void* p); //For shell integration.
		static double btb_wildlife(double d, double mean, double stdev); //Normal distribution.
		static double btb_wildlife_shell(double d, void* p); //For shell integration.
		void normalize_btb_kernels();


	public:
		///> Constructs a kernel from an equation (determined by variable kType)
		Local_spread(int kernelType,
			std::vector<std::vector<double>> kparams = std::vector<std::vector<double>>());
		///> Constructs a kernel from external file
		Local_spread(int kernelType,
			std::string fname);
		~Local_spread();
		///> Calculates or matches kernel value according to form defined at construction
		double atDistSq(double);
		///> Calculates kernel value according to form with multiple arguments defined at construction. Used for bTB.
		double atDistSq(double, double, int);
};

#endif // Local_spread_h
