#ifndef Control_resource_h
#define Control_resource_h

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include "Farm.h"

/// Struct containing location information for a Control_resource
struct resourceLocation
{
	double x;
	double y;
	std::string fips;
	std::string stateAbbrev;
};


/// Object for a type of control resource (i.e. single landfill or all landfills in a state).
class Control_resource
{
	protected:				
		double x_coordinate;
		double y_coordinate;
		std::string fips;
		std::string state;
//  cell id could also be stored here
			
		double capacity; ///< Estimate of number of animals that can be treated using this resource
		std::tuple<double, double> fixedDailyLimit;
//	std::unordered_map<int, int> dailyLimits; //Can specify per-day capacity changes starting at time (index)
		
	public:
		Control_resource(resourceLocation, int startLevel_in=0);
		Control_resource();
		~Control_resource();	
		
		void set_fixedDailyLimit(std::tuple<double, double>); //inlined
		int get_dailyLimit() const;
		double get_capacity() const;
		void add_capacity(double); //inlined
		void subtract_capacity(double toSubtract = 1); //inlined
};
/// Set value (during initial file read-in) to be used as fixed daily limit
inline void Control_resource::set_fixedDailyLimit(std::tuple<double, double> limit)
{
	fixedDailyLimit = limit;
}
///
inline double Control_resource::get_capacity() const
{
	return capacity;
}
/// Add capacity (i.e. for loading a particular landfill as part of a state Control_resource)
inline void Control_resource::add_capacity(double toAdd)
{
	capacity += toAdd;
}

/// Subtract capacity (as resource capacity diminishes)
inline void Control_resource::subtract_capacity(double toSubtract/*=1*/)
{
	capacity -= toSubtract;
}

#endif //Control_resource_h