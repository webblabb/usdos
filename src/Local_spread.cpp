#include "Local_spread.h"

///	\param[in]	kernelType	Integer specifying type of equation-based kernel. 0 uses 
///	\f$\frac{k_1}{1+(\frac{d}{k_2})^{k_3}}\f$. 2 uses \f$\frac{k_1}{(1+\frac{d}{k_2})^{k_3}}\f$.
Local_spread::Local_spread(int kernelType, std::vector<double> kparams)
	:
	kType(kernelType),
	kp(kparams)
{	
	verbose = verboseLevel;
	switch (kType)
	{
		case 0:{ // power law function
		// already set: kp[0] is k1, [1] is k2, [2] is k3
		// precalculate k3/2 and add as fourth parameter at [3]
		// precalculate k2^k3 and add as fifth parameter at [4]
			kp.emplace_back(kp.at(2)/2);
			kp.emplace_back(pow(kp.at(1),kp.at(2)));
			break;
		}
		case 2:{ // "kernel4" from JapanFMD
		// use parameters as entered
			break;
		}
		default:{
			std::cout << "Kernel type must be 0 or 2 if constructing with parameters. Exiting..." << std::endl; exit(EXIT_FAILURE);
		}
	}
}

Local_spread::Local_spread(int kernelType, std::string fname)
	:
	kType(kernelType),
	datafile(fname)
{	
	verbose = verboseLevel;
	switch (kType)
	{
		case 1:{ // data-based levels
std::cout<<"Creating a kernel with datafile "<<datafile<<std::endl;
		// read in levels from file
			std::ifstream d(datafile);
			if(!d){std::cout << "Data-based local spread file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
			if(d.is_open()){
if(verbose>1){std::cout << "Data-based local spread file open." << std::endl;}
				while(! d.eof()){
					std::string line;
					getline(d, line); // get line from file "f", save as "line"
					std::vector<std::string> line_vector = split(line, ' '); // separate by tab
					if(! line_vector.empty()){ // if line_vector has something in it
						double distM = stringToNum<double>(line_vector[0]);
						// store distance as distance-squared
						distProb[distM*distM] = stringToNum<double>(line_vector[1]);
if(verbose>1){std::cout<<" At "<<distM<<", risk is "<<distProb[distM*distM]<<std::endl;}
					} // close "if line_vector not empty"
				} // close "while not end of file"
			} // close "if file is open"
			break;
		}
		default:{
			std::cout << "Kernel type must be 1 if constructing from file. Exiting..." << std::endl; exit(EXIT_FAILURE);
		}
	}
}

Local_spread::~Local_spread()
{
}

/// \param[in] distSq	distance squared
/// \returns kernel value at distance
///	 
/// For kType 0: to illustrate that the form of kernel(dist squared):
/// \f$\frac{k_1}{1+\frac{dsq^{k_3/2}}{k_2^{k_3}}}\f$
/// gives the same output as kernel(dist):
/// \f$\frac{k_1}{1+\frac{d}{k_2}^{k_3}}\f$, 
/// run in R:
/// \code{.py}
/// usedist = 1:2000
/// k1 = 0.089
///	k2 = 1000
/// k3 = 3
/// # original f1(distance)
/// plot(k1 / (1 + (usedist/k2)^k3) ~ usedist) 
/// usq = usedist^2
/// # add points for f2(distance-squared)
/// points(usedist, (k1 / (1 + (usq^(k3/2))/(k2^k3))),col="blue",pch="*")
/// \endcode

double Local_spread::atDistSq(double distSq)
{
	double k;
	switch (kType)
	{
		case 0:{ // power law function
			k = kp.at(0)/(1+pow(distSq,kp.at(3))/kp.at(4));
			break;
		}
		case 1:{ // UK data-based levels
			// if distance is past last value, use last value
			if (distSq >= distProb.rbegin()->first){k = 0;			
if(verbose>1){std::cout << "Distance^2 past end, set to 0"<<std::endl;}
			} else {
				auto equalOrHigher = distProb.lower_bound(distSq); // finds equal or higher match to distSq key
				k = equalOrHigher->second;
				// check if the previous key (if applicable) is closer
				if (equalOrHigher != distProb.begin()){
					auto lower = std::prev(equalOrHigher); // gets one key before that
					// use whichever has less difference (is closer)
					if ( std::abs(distSq - lower->first) < std::abs(equalOrHigher->first - distSq)){
						k = lower->second;
					}
				}
if(verbose>1){std::cout << "Matched distance "<<std::sqrt(distSq)<<" with probability value of "<<k<<std::endl;}
			}
			break;
		}
		case 2:{ // "kernel4"
			double dOverK2 = std::sqrt(distSq)/kp.at(1);
			k = kp.at(0)/pow((1+dOverK2),kp.at(2));
			break;
		}
		default:{
			std::cout << "Unrecognized kernel type. Exiting..." << std::endl; exit(EXIT_FAILURE);
		}
	}
 return std::min(1.0,k);
}