//© 2024 Colorado State University
#include <Shipment_kernel.h>
#include <County.h>
#include <shared_functions.h>
#include <iostream>

//Variable to store binned distances in, static since it is identical to all shipment_kernel objects.
std::vector<double> Shipment_kernel::binned_distances;

Shipment_kernel::Shipment_kernel(double a, double b, std::string type, bool binning_on) :
    a(a), b(b), binning_on(binning_on)
{
    if(type.compare("USAMMv1") == 0)
    {
        k_function = &Shipment_kernel::linear_distance_kernel;
        d_function = &Shipment_kernel::linear_euclidean;
        binning_function = &Shipment_kernel::get_bin_USAMMv2;
    }
    else if(type.compare("USAMMv3") == 0)
    {
        k_function = &Shipment_kernel::local_kernel_from_half_and_R_func;
        d_function = &Shipment_kernel::linear_euclidean;
        binning_function = &Shipment_kernel::get_bin_USAMMv3;
    }
    else if(type.compare("quadratic") == 0)
    {
        k_function = &Shipment_kernel::quadratic_distance_kernel;
        d_function = &Shipment_kernel::quadratic_euclidean;
        binning_function = &Shipment_kernel::get_bin_USAMMv2;
    }
    else if(type.compare("power_exp") == 0)
    {
        k_function = &Shipment_kernel::power_exp_from_half_and_R_func;
        d_function = &Shipment_kernel::linear_euclidean;
        binning_function = &Shipment_kernel::get_bin_USAMMv2;
    }
    else if(type.compare("one_minus_exp") == 0)
    {
        k_function = &Shipment_kernel::one_minus_exp_from_half_and_R_func;
        d_function = &Shipment_kernel::linear_euclidean;
        binning_function = &Shipment_kernel::get_bin_USAMMv2;
    }
    else if(type.compare("local") == 0)
    {
        k_function = &Shipment_kernel::local_kernel_from_half_and_R_func;
        d_function = &Shipment_kernel::linear_euclidean;
        binning_function = &Shipment_kernel::get_bin_USAMMv2;
    }
    else
    {
        std::cout <<"Error: Unknown shipment kernel type: " << type << ". Exiting..."
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    a_sq = a*a;
    b_half = b * 0.5;

    if(binning_on and binned_distances.empty())
    {
        set_bins_peters(50);
    }
}

Shipment_kernel::~Shipment_kernel() {}

double Shipment_kernel::kernel(County* c1, County* c2)
{
    double d = (this->*d_function)(c1, c2);
    return (this->*k_function)(d);
}

double Shipment_kernel::distance(County* c1, County* c2, std::string type)
{
    if(type == "linear")
    {
        return linear_euclidean(c1, c2);
    }

    else if(type == "quadratic")
    {
        return quadratic_euclidean(c1, c2);
    }
    else
    {
        return (this->*d_function)(c1, c2);
    }
}

void Shipment_kernel::set_bin_size(double in_size)
{
    bin_size = in_size;
}

void Shipment_kernel::set_longest_distance(double in_dist)
{
    longest_distance = in_dist;
    if(binning_on)
    {
        set_bins_peters(50); //Update bins with the new distance.
    }
}

void Shipment_kernel::set_bins_unif()
{
    binned_distances.clear();
    int n_bins = int(longest_distance / bin_size) + 1;
    binned_distances.reserve(n_bins);
    for(int i = 0; i < n_bins; i++)
    {
        binned_distances.push_back(i*bin_size);
    }
}

void Shipment_kernel::set_bins_peters(int no_sub_ints){

    double interval_len = longest_distance / double(no_sub_ints);
    std::vector<double> sub_vec;
    std::vector<double> bin_limits;

    for(int i = 0; i < no_sub_ints-1; ++i)
    {
        sub_vec.push_back((i+1) * interval_len);
    }
    sub_vec.push_back(longest_distance);

    // Add more intervals in the beginning of the vector
    for(int i=1;i<(no_sub_ints);i++)
    {
        bin_limits.push_back((i)*interval_len / ((no_sub_ints)));
    }

    for(int j = 0; j < static_cast<int>(no_sub_ints / 2.0); j++)
    {
        for(int i = 1; i < (no_sub_ints-j); i++)
        {
            bin_limits.push_back(sub_vec[j]+(i)*interval_len / ((no_sub_ints-j)));

        }
    }

    for(int j = static_cast<int>(no_sub_ints / 2); j < no_sub_ints - 1; j++)
    {
        for(int i = 0; i < static_cast<int> (no_sub_ints / 2); i++)
        {
            bin_limits.push_back(sub_vec[j] + (i)*interval_len / static_cast<int>(no_sub_ints / 2));
        }
    }

    bin_limits.push_back(longest_distance + 1);
    bin_limits.push_back(longest_distance + 2);
    size_t no_dist_ints = bin_limits.size();
    binned_distances.resize(no_dist_ints);
    binned_distances[0] = bin_limits[0] / 2.0;

    for(size_t i = 1; i < no_dist_ints; i++)
    {
        binned_distances[i] = bin_limits[i-1] + (bin_limits[i] - bin_limits[i-1]) / 2.0;
    }
}

double Shipment_kernel::get_bin_USAMMv2(double d)
{
    //Binary search for the correct bin of d
    bool done = false;
    int lower = 0; //Lower index of current sub-vector
    int upper = binned_distances.size() - 1; //Upper index of current sub-vector
    int mid = 0;
    int n_tries = 0;

    do
    {
        mid = lower + int(ceil((upper-lower) / 2)); //Get mid index of current sub-vector.
        if(d == binned_distances[mid])
        {
            done = true;
            return binned_distances[mid];
        }
        else if(d > binned_distances[mid])
        {
            lower = mid;
        }
        else if(d < binned_distances[mid])
        {
            upper = mid;
        }
        if(abs(upper - lower) == 1)
        {
            double udiff = std::abs(binned_distances[upper] - d);
            double ldiff = std::abs(binned_distances[lower] - d);
            if(ldiff < udiff)
            {
                done = true;
                return binned_distances[lower];
            }
            else
            {
                done = true;
                return binned_distances[upper];
            }
        }
        n_tries += 1;
        if(n_tries >= 100)
        {
            std::cout << "Stuck in when getting bin for distance " << d << ". Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
    } while(!done);

    return -1;
}

double Shipment_kernel::get_bin_USAMMv3(double d)
{
    d /= 1000.0; //To km.
    if(d < 0 )
    {
        std::cout << "Negative distance in distance binning function." << std::endl;
        exit(EXIT_FAILURE);
    }
    else if(d >= 0 and d < 100)
    {
        return 1000.0 * int(d + 0.5);
    }
    else if(d >= 100 and d < 1000)
    {
        return 1000.0 * int((d / 10) + 0.5) * 10;
    }
    else
    {
        // d >= 1000
        return 1000.0 * int((d / 100) + 0.5) * 100;
    }
}

double Shipment_kernel::linear_distance_kernel(double d)
{
    return std::exp(-std::pow(d/a,b));
}

double Shipment_kernel::quadratic_distance_kernel(double sq_d)
{
    return std::exp(-std::pow(sq_d / a_sq, b_half));
}

//These comments are relevant for the three kernel implementations
//one_minus_exp_from_half_and_R_func, power_exp_from_half_and_R_func (same
//as original kernel) and local_kernel_from_half_and_R_func:
// d_half=40;%distance where function has dropped to x1
// R=50;% Ratio between distance where function is at x1 and distance where it is at x2.
// x1=.5; %First relevant function value (typically 0.5)
// x2=.05;%Second relevant function value (typically 0.05 or 0.01)
double Shipment_kernel::power_exp_from_half_and_R_func(double d)
{
  return std::exp(-std::pow(d/a, b));
}

double Shipment_kernel::one_minus_exp_from_half_and_R_func(double d)
{
    return oneMinusExp(-std::pow(d/a, b));
}

double Shipment_kernel::local_kernel_from_half_and_R_func(double d)
{
    return 1 / (1 + std::pow(d/a, b));
}

double Shipment_kernel::linear_euclidean(County* c1, County* c2)
{
    //If origin and destination is the same county the distance is set to the
    //average distance between two random points within a square with the
    //same area as the county (which is side of square * 0.5214).
    double d;
    if(c1 == c2)
    {
        d = std::sqrt(c1->get_area()) * 0.5214;
    }
    else
    {
        const Point* p1 = c1->get_centroid();
        const Point* p2 = c2->get_centroid();
        d = sqrt((p1->x - p2->x) * (p1->x - p2->x) +
                 (p1->y - p2->y) * (p1->y - p2->y));
    }

    if(binning_on)
    {
        d = (this->*binning_function)(d);
    }
    return d;
}

double Shipment_kernel::quadratic_euclidean(County* c1, County* c2)
{
    double d;
    if(c1 == c2)
    {
        d = c1->get_area() * 0.5214 * 0.5214;
    }
    else
    {
        const Point* p1 = c1->get_centroid();
        const Point* p2 = c2->get_centroid();
        d = ((p1->x - p2->x) * (p1->x - p2->x) +
            (p1->y - p2->y) * (p1->y - p2->y));
    }

    if(binning_on)
    {
        d = (this->*binning_function)(d);
    }
    return d;
}
