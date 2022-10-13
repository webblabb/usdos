#ifndef KERNEL_F_H
#define KERNEL_F_H

#include <vector>
#include <string>
#include <iostream>

/*!
A kernel function for shipments. Construct with parameters a and b
as well as a string describing the type, either "linear" or
"quadratic". The a and b parameters control the shape of the kernel
function.
"linear" uses a standard kernel that is a function of
distance, d. The "quadratic" alternative uses a kernel that is a
function of d^2 to avoid square root operations. Both give
identical results. To get a kernel value for a pair of counties
call the function kernel(c1, c2), where the arguments are pointers to
two counties. If both pointers point to the same county object
the kernel value that is returned is for d = the average distance
within a square of the same area as the county. The functions for
getting the distance between two counties are built into the class.
Also, for convenience, the function distance(c1, c2, type) can
be used to get the distance between two counties. The type argument
defaults to "linear" so if the actual euclidean distance is desired
the type argument can be omitted.
The shipment kernel class also manages the binning of distances. The
binning feature is optional and does not improve performance. The reason
for binning is that USAMM requires binning of distances. This means that
the parameters used for shipping generation are estimated using binned
distances and so the distances used for kernel evaluation in the
simulation should be binned in the same way.
*/

class County;

class Shipment_kernel
{
typedef double (Shipment_kernel::*k_fun_ptr)(double); //Kernel function pointer
typedef double (Shipment_kernel::*d_fun_ptr)(County*, County*); //Distance function pointer
typedef double (Shipment_kernel::*bin_fun_ptr)(double); //Binning function pointer.

public:
    Shipment_kernel(double dx, double R, std::string type = "linear", bool binning_on = false);
    ~Shipment_kernel();

    double kernel(County* c1, County* c2);
    ///Calls the correct distance function (depending on what type this distance kernel object is of)
    ///and returns the distance between two county objects.
    double distance(County* c1, County* c2, std::string type = "linear");
    ///Sets the size of each bin
    void set_bin_size(double in_size);
    ///Sets the maximum distance to bin for.
    void set_longest_distance(double in_dist);

private:
    double a, b;
    bool binning_on;
    double a_sq, b_half;
    double R, dx1; //These parameters are used to get the points at which the kernel value = x1 (dx1) and value at x1 / value at x2.
    double bin_size = 20000;
    double longest_distance = 6000000;
    static std::vector<double> binned_distances;
    k_fun_ptr k_function;
    d_fun_ptr d_function;
    bin_fun_ptr binning_function;

    ///Set bins to be smaller at short distances and increase gradually in size
    ///with the distance. This is the way that binning is implemented in USAMM.
    void set_bins_peters(int no_sub_ints);
    ///Set bins to be uniformly spaced.
    void set_bins_unif();
    ///Finds the bin of the distance d given the current set of bins.
    double get_bin_USAMMv2(double d);
    ///Takes a distance d and returns it binned as needed for USAMMv3 shipments.
    double get_bin_USAMMv3(double d);
    ///The distance kernel function of USAMM.
    double linear_distance_kernel(double d);
    ///The distance kernel function of the squared distance. Experimental. Possibly faster.
    double quadratic_distance_kernel(double sq_d);
    ///Additional kernel 1
    double one_minus_exp_from_half_and_R_func(double d);
    ///Same as 'original' kernel (linear_distance_kernel())
    double power_exp_from_half_and_R_func(double d);
    ///Additional kernel 2
    double local_kernel_from_half_and_R_func(double d);
    ///Normal euclidean distance between two counties.
    double linear_euclidean(County* c1, County* c2);
    ///Returns the squared distance between two counties. See class
    ///description for more info.
    double quadratic_euclidean(County* c1, County* c2);

};
#endif // KERNEL_F_H
