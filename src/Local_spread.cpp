#include "Local_spread.h"
#include "gsl/gsl_integration.h"

///	\param[in]	kernelType	Integer specifying type of equation-based kernel. 0 uses
///	\f$\frac{k_1}{1+(\frac{d}{k_2})^{k_3}}\f$. 2 uses \f$\frac{k_1}{(1+\frac{d}{k_2})^{k_3}}\f$.
Local_spread::Local_spread(int kernelType, std::vector<std::vector<double>> kparams)
	:
	kType(kernelType)

{
    kp_q1 = kparams.at(0);
    if(kparams.size() == 4)
    {
        kp_q2 = kparams.at(1);
        kp_q3 = kparams.at(2);
        kp_q4 = kparams.at(3);
    }
    else if(kparams.size() != 1)
    {
        std::cout << "The local spread kernel parameters must either be yearly (one set of parameters) or quarterly (four sets)." << std::endl;
        exit(EXIT_FAILURE);
    }


	verbose = verboseLevel;
	switch (kType)
	{
		case 0:{ // power law function
		// already set: kp_q1[0] is k1, [1] is k2, [2] is k3
		// precalculate k3/2 and add as fourth parameter at [3]
		// precalculate k2^k3 and add as fifth parameter at [4]
			kp_q1.emplace_back(kp_q1.at(2)/2);
			kp_q1.emplace_back(pow(kp_q1.at(1),kp_q1.at(2)));
			break;
		}
		case 2:{ // "kernel4" from JapanFMD
		// use parameters as entered
			break;
		}
		case 3:{ //bTb kernel.
		    // kparams by idx are 0: phi (local), 1: alpha (local), 2: psi (wildlife), 3: beta (wildlife), 4: x (relative weight of wildlife [0,1]), 5: scaling factor for entire combined kernel applied after normalization.
            normalize_btb_kernels();
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
			k = kp_q1.at(0)/(1+pow(distSq, kp_q1.at(3))/kp_q1.at(4));
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
			double dOverK2 = std::sqrt(distSq)/kp_q1.at(1);
			k = kp_q1.at(0)/pow((1+dOverK2),kp_q1.at(2));
			break;
		}
		default:{
			std::cout << "Unrecognized kernel type. Exiting..." << std::endl; exit(EXIT_FAILURE);
		}
	}
 return std::min(1.0,k);
}

//Additional path for kernels that have three independent dimensions (e.g. distance, wildlife density, quarter)
double Local_spread::atDistSq(double distSq, double kernel_dim_two, int quarter_idx)
{
	double k;
	switch (kType)
	{
		case 3:{
		    // kparams by idx are 0: phi (local), 1: alpha (local), 2: psi (wildlife), 3: beta (wildlife), 4: x (relative weight of wildlife [0,1]), 5: scaling factor for entire combined kernel applied after normalization.
            double d = std::sqrt(distSq);
            switch (quarter_idx)
            {
                case 0:
                    k = kp_q1[5] * btb_mixed_kernel(d, kernel_dim_two, kp_q1[0], kp_q1[1], kp_q1[2], kp_q1[3], kp_q1[4], quarter_idx);
                    break;
                case 1:
                    k = kp_q2[5] * btb_mixed_kernel(d, kernel_dim_two, kp_q2[0], kp_q2[1], kp_q2[2], kp_q2[3], kp_q2[4], quarter_idx);
                    break;
                case 2:
                    k = kp_q3[5] * btb_mixed_kernel(d, kernel_dim_two, kp_q3[0], kp_q3[1], kp_q3[2], kp_q3[3], kp_q3[4], quarter_idx);
                    break;
                case 3:
                    k = kp_q4[5] * btb_mixed_kernel(d, kernel_dim_two, kp_q4[0], kp_q4[1], kp_q4[2], kp_q4[3], kp_q4[4], quarter_idx);
                    break;
                default:{
                    std::cout << "Unrecognized kernel type. Exiting..." << std::endl; exit(EXIT_FAILURE);
                }
            }
            break;
		}
		default:{
			std::cout << "Unrecognized kernel type. Exiting..." << std::endl; exit(EXIT_FAILURE);
		}
	}
	if(std::isnan(k) or std::isinf(k))
    {
        std::cout << "qwdqwd" << std::endl;
    }
    return k;
}

double Local_spread::btb_mixed_kernel(double d, double wildl_dens, double lp1,
                                      double lp2, double wp1, double wp2, double x,
                                      int quarter_idx)
{
    //Mixture kernel consisting of local (L) and wildlife (W) transmission.
    //Both kernels are normalized by their volume over d = [0, 1.0e6] m, so that
    //the volume of both are == 1.0.
    //The proportional influence of the wildlife kernel is controlled by x. When
    //x = .5, both kernels have equal weight.
    //After kernels have been evaluated, the combined effect is scaled up so that
    //the volume of the entire mixed kernel is equal to the volume of the local
    //kernel in order to get comparable infection pressure to VanderWaal et al.
    //lp1, lp2 are local kernel parameters 1 & 2; wp1, wp2 are wildl. kernel parameters 1 & 2.
    return (x * wildl_dens * btb_wildlife(d, wp1, wp2) / V_btb_wildlife[quarter_idx] +
           (1-x) * btb_local(d, lp1, lp2) / V_btb_local[quarter_idx]) * V_btb_local[quarter_idx];
}

double Local_spread::btb_local(double d, double lp1, double lp2)
{
    return lp1 * std::exp(-d*lp2);
}

double Local_spread::btb_local_shell(double d, void* p)
{
    std::vector<double>* btb_pars = (std::vector<double>*)p;
    double phi = (*btb_pars)[0];
    double alpha = (*btb_pars)[1];
    return d * btb_local(d, phi, alpha);
}

double Local_spread::btb_wildlife(double d, double mean, double stdev)
{
    return gsl_ran_gaussian_pdf(d-mean, stdev);
}

double Local_spread::btb_wildlife_shell(double d, void* p)
{
    std::vector<double>* btb_pars = (std::vector<double>*)p;
    double mean = (*btb_pars)[2];
    double stdev = (*btb_pars)[3];
    return d * btb_wildlife(d, mean, stdev);
}

void Local_spread::normalize_btb_kernels()
{
    // kparams by idx are 0: phi (local), 1: alpha (local), 2: psi (wildlife), 3: beta (wildlife), 4: x (relative weight of wildlife [0,1]), 5: scaling factor for entire combined kernel applied after normalization.
    //Begin with calculating scaling factors to scale local and wildlife components so they have volume=1.0 respectively.
    //We'll integrate over the interval 0,1000 km as (K(d) = 1.0e-324 at around d = 508 km with the parameters given
    //in VanderWaal et al. 2017.

    std::vector<std::vector<double>*> ptr_vec = { &kp_q1, &kp_q2, &kp_q3, &kp_q4 };
    V_btb_local.resize(ptr_vec.size(), -1.0);
    V_btb_wildlife.resize(ptr_vec.size(), -1.0);
    for(size_t i=0; i<ptr_vec.size(); ++i)
    {

        size_t ws_size = 1000;
        gsl_integration_workspace * ws = gsl_integration_workspace_alloc(ws_size);

        double result_local, error_local;
        gsl_function F_local;
        F_local.function = &Local_spread::btb_local_shell;
        F_local.params = ptr_vec[i];
        gsl_integration_qag(&F_local, 0.0, 1000.0*1000.0, 1e-7, 1000.0,
                            ws_size, 6, ws, &result_local, &error_local);

        double result_wildlife, error_wildlife;
        gsl_function F_wildlife;
        F_wildlife.function = &Local_spread::btb_wildlife_shell;
        F_wildlife.params = ptr_vec[i];
        gsl_integration_qag(&F_wildlife, 0.0, 1000.0*1000.0, 1e-7, 1000.0,
                            ws_size, 6, ws, &result_wildlife, &error_wildlife);

        //After integrating the kernels over the interval [0, 1.0e6] meters, we calculate the volumes:
        double local_volume = result_local * 2 * M_PI; //Volume by shell integration (rotation around y).
        if(local_volume <= 0.0 or std::isnan(local_volume) or std::isinf(local_volume))
        {
            std::cout << "Failed to find volume under the local btb kernel for quarter " << i+1
                      << "; volume = " << local_volume << std::endl;
            exit(EXIT_FAILURE);
        }
        V_btb_local[i] = local_volume;

        double wildl_volume = result_wildlife * 2 * M_PI;
        if(wildl_volume <= 0.0 or std::isnan(wildl_volume) or std::isinf(wildl_volume))
        {
            std::cout << "Failed to find volume under the btb wildlife kernel for quarter " << i+1
                      << "; volume = " << wildl_volume << std::endl;
            exit(EXIT_FAILURE);
        }
        V_btb_wildlife[i] = wildl_volume;

        gsl_integration_workspace_free (ws);
    }
}
