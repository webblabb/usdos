//© 2024 Colorado State University
#include <cmath>
#include <iostream>
#include <algorithm>
#include <Region.h>

Region::Region(std::string id) :
    id(id), centroid(0.0, 0.0)
{
    set_initialized(is_set_id);
}

Region::Region(std::string id, double x, double y) :
    id(id), centroid(x, y)
{
    set_initialized(is_set_id);
    set_initialized(is_set_position);
}

Region::~Region() {}

double Region::distance_to(Region* target) const
{
    return measure_distance(&centroid, target->get_centroid());
}

double Region::measure_distance(const Point* p1, const Point* p2) const
{
    return std::sqrt(std::pow(p1->x - p2->x, 2) +
                     std::pow(p1->y - p2->y, 2));
}

void Region::set_position(double x, double y)
{
    centroid.x = x;
    centroid.y = y;
    set_initialized(is_set_position);
}

void Region::set_initialized(bool& parameter)
{
    parameter = true;
    all_initialized();
}

void Region::all_initialized()
{
    if(is_set_id and is_set_position)
    {
        region_initialized = true;
    }
}

void Region::not_initialized()
{
    std::cout << "Error: " << type << " " << id << " has not yet been completely initialized." << std::endl
              << "Exiting...";
    exit(EXIT_FAILURE);
}

Region_status::Region_status(Region* r)
	:
	Region(*r), // call copy constructor of Region to fill in all other members
	reported(false)
{
	probPreventExposure.emplace_back(0);
}

Region_status::~Region_status()
{
}

int Region_status::get_start(std::string s) const
{
	int output = -1;
	if (start.count(s) > 0){
		output = start.at(s);
	}
	return output;
}

int Region_status::get_end(std::string s) const
{
	int output = -1;
	if (end.count(s) > 0){
		output = end.at(s);
	}
	return output;
}

/// Returns the probability that exposure is prevented, based on the
/// effectiveness of currently effective control types. If multiple control types are
/// currently effective, the maximum effectiveness value is returned (no additive
/// effects are assumed).
double Region_status::get_probPreventExposure() const
{
	auto ppe = std::max_element(probPreventExposure.begin(), probPreventExposure.end());
	double output = *ppe;
	if (output > 1){output = 1;
	} else if (output < 0){output = 0;}
	return output;
}

void Region_status::rem_probPreventExposure(double eff)
{
	auto it = std::find(probPreventExposure.begin(), probPreventExposure.end(), eff);
	if (it != probPreventExposure.end()) { // if value is found
		// swap the one to be removed with the last element
		std::vector<double>::iterator last = std::prev(probPreventExposure.end());
		std::iter_swap(it, last);
		// and remove the item at the back of the container
		probPreventExposure.pop_back();
	}
}

/// Returns the probability that transmission is prevented, based on the
/// effectiveness of currently effective control types. If multiple control types are
/// currently effective, the maximum effectiveness value is returned (no additive
/// effects are assumed).
double Region_status::get_probPreventTransmission() const
{
	auto ppe = std::max_element(probPreventTransmission.begin(), probPreventTransmission.end());
	double output = *ppe;
	if (output > 1){output = 1;
	} else if (output < 0){output = 0;}
	return output;
}

void Region_status::rem_probPreventTransmission(double eff)
{
	auto it = std::find(probPreventTransmission.begin(), probPreventTransmission.end(), eff);
	if (it != probPreventTransmission.end()) { // if value is found
		// swap the one to be removed with the last element
		std::vector<double>::iterator last = std::prev(probPreventTransmission.end());
		std::iter_swap(it, last);
		// and remove the item at the back of the container
		probPreventTransmission.pop_back();
	}
}
