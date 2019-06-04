#include <iostream>
#include <vector>
#include <algorithm> // for std::max_element (susceptibility/infectiousness)

#include "Grid_cell.h"

/// Constructed with cell dimensions and premises within. Calculates and stores maximum
/// transmission values (for overestimating transmission probabilities).
Grid_cell::Grid_cell(const int in_id, const double in_x, const double in_y,
	const double in_s, const std::vector<Farm*> in_farms)
	:
	id(in_id),
	x(in_x),
	y(in_y),
	s(in_s),
	farms(in_farms)
{
	// Add all farms' susceptibility and infectiousness to respective vectors and find max
	std::vector <double> allSus;
	std::vector <double> allInf;
	for (auto& f:farms){
		allSus.emplace_back(f->Farm::get_sus_max()); // add farm's susceptibility to vector
		allInf.emplace_back(f->Farm::get_inf_max()); // add farm's infectiousness to vector
		// Could make this optional if control or specific target type turned off
		// Add counties & states
		countiesIncluded.emplace(f->Farm::get_parent_county()->Region::get_id());
		statesIncluded.emplace(f->Farm::get_parent_county()->Region::get_id());
		}
	maxSus = *std::max_element(allSus.begin(),allSus.end());
	maxInf = *std::max_element(allInf.begin(),allInf.end());

	// Calculate boundaries of each side
	south = y; // lower boundary
	north = y+s; // upper boundary
	west = x; // leftmost boundary
	east = x+s; // rightmost boundary

}

Grid_cell::~Grid_cell()
{
}

void Grid_cell::addNeighbor(Grid_cell* in_neighbor)
{
	neighbors.emplace_back(in_neighbor);
}

///> Swaps contents of the empty member map for calculated kernel*susceptibility values
void Grid_cell::take_KernelValues(std::unordered_map<int, double>& in_kern)
{
	susxKern.swap(in_kern);
}

void Grid_cell::removeFarmSubset(std::vector<int>& toRemove)
{
	auto newEnd = std::remove_if(farms.begin(),farms.end(),farmIDpresent(toRemove));
	std::vector<Farm*> newFarms(farms.begin(),newEnd);
	farms.swap(newFarms);
}
