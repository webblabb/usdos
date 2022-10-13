#ifndef Diagnostic_resource_h
#define Diagnostic_resource_h

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include "Farm.h"



/// Object for a type of diagnostic resource (i.e. number of PCR kits).
class Diagnostic_resource
{
	protected:				
		int diagnosticCapacity; ///< Estimate of number of animals that can be tested using this resource
		std::tuple<double, double> diagnosticFixedDailyLimit;
//	std::unordered_map<int, int> dailyLimits; //Can specify per-day capacity changes starting at time (index)
		
	public:
		Diagnostic_resource(int startLevel_in=0);
		Diagnostic_resource();
		~Diagnostic_resource();	
		
		void set_diagnosticFixedDailyLimit(std::tuple<double, double>); //inlined
		int get_diagnosticDailyLimit() const;
		int get_diagnosticCapacity() const;
		void add_diagnosticCapacity(int); //inlined
		void subtract_diagnosticCapacity(int toSubtract = 1); //inlined
};
/// Set value (during initial file read-in) to be used as fixed daily limit
inline void Diagnostic_resource::set_diagnosticFixedDailyLimit(std::tuple<double, double> limit)
{
	diagnosticFixedDailyLimit = limit;
}
///
inline int Diagnostic_resource::get_diagnosticCapacity() const
{
	return diagnosticCapacity;
}
/// Add capacity 
inline void Diagnostic_resource::add_diagnosticCapacity(int toAdd)
{
	diagnosticCapacity += toAdd;
}
/// Subtract capacity (as resource capacity diminishes)
inline void Diagnostic_resource::subtract_diagnosticCapacity(int toSubtract/*=1*/)
{
	diagnosticCapacity -= toSubtract;
}

#endif //Diagnostic_resource_h