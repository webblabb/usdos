#ifndef grid_cell_h
#define grid_cell_h

#include "County.h" // to include counties included in each cell
#include "shared_functions.h" // for isWithin function in struct farmInList
#include "Farm.h"
#include "State.h"

/// Each Grid_cell has x-y coordinates of lower left corner, max susceptibility/infectiousness,
/// vectors of premises within, neighbors, kernel values to other cells according to
/// its most infectious premises, boundaries to the north, south, east, west, and states
/// included in the cell.
class Grid_cell
{
    private:
    	int id; /// Unique identifier for this cell
        double x; /// x-coordinate of lower left corner of cell (same units as Farms)
        double y; /// y-coordinate of lower left corner of cell (same units as Farms)
        double s; /// length of one side of square cell (same units as x,y)
        double maxSus; /// Maximum susceptibility value of all premises in this cell
        double maxInf; /// Maximum infectiousness value of all premises in this cell
        std::vector<Farm*> farms; /// Maximum susceptibility value of all premises in this cell
        std::vector<Grid_cell*> neighbors; /// All Grid_cells touching this cell, not including self
        std::unordered_map<int, double> susxKern; /// Map of pre-calculated values for cell-cell maximum susceptibility * distance-based kernel. Map key is int rather than pointer because cells are copied and modified in Grid_checker.
				std::set<std::string> statesIncluded; /// States (2 letter abbreviation) included in cell
        std::set<std::string> countiesIncluded; /// Counties (identified by FIPS code as string) included in this cell

    		double north; /// y-coordinate of upper corners of cell (same units as Farm coordinates)
    		double south; /// y-coordinate of lower corners of cell (same units as Farm coordinates)
    		double east; /// x-coordinate of right corners of cell (same units as Farm coordinates)
    		double west; /// x-coordinate of left corners of cell (same units as Farm coordinates)
    		
    public:
		Grid_cell(const int, const double, const double, const double, const std::vector<Farm*>);
		~Grid_cell();

        void addNeighbor(Grid_cell*);
        std::vector<Farm*> get_farms() const; //inlined
		int get_id() const; //inlined
        double get_maxInf() const; //inlined
        double get_maxSus() const; //inlined
		const std::vector<Grid_cell*>* get_neighbors(); //inlined
        double get_num_farms() const; // inlined
        double get_s() const; // inlined
        const std::unordered_map<int, double>* get_susxKernel(); //inlined
        double get_x() const; // inlined
        double get_y() const; // inlined
        double get_south() const; // inlined
        double get_north() const; // inlined
        double get_east() const; // inlined
        double get_west() const; // inlined
        std::set<std::string> get_counties() const; //inlined
        std::set<std::string> get_states() const; //inlined
        double kernelTo(int) const; //inlined
        void removeFarmSubset(std::vector<int>&);
		void take_KernelValues(std::unordered_map<int, double>&);

};

inline std::vector<Farm*> Grid_cell::get_farms() const {
	return farms;}

inline int Grid_cell::get_id() const {
	return id;}

inline double Grid_cell::get_maxInf() const {
	return maxInf;}

inline double Grid_cell::get_maxSus() const {
	return maxSus;}

inline const std::vector<Grid_cell*>* Grid_cell::get_neighbors(){
	return &neighbors;}

inline double Grid_cell::get_num_farms() const {
	return double(farms.size());}

inline double Grid_cell::get_s() const {
    return s;}

inline double Grid_cell::get_x() const {
	return x;}

inline double Grid_cell::get_y() const {
    return y;}

inline double Grid_cell::get_south() const {
    return south;}

inline double Grid_cell::get_north() const {
    return north;}
    
inline double Grid_cell::get_east() const {
    return east;}
    
inline double Grid_cell::get_west() const {
    return west;}

inline std::set<std::string> Grid_cell::get_counties() const {
		return countiesIncluded;}
		
inline std::set<std::string> Grid_cell::get_states() const {
		return statesIncluded;}

inline double Grid_cell::kernelTo(int id) const {
	return susxKern.at(id);}

///> Identifies if a Farm* is present in a list of ID numbers
struct farmIDpresent // used in removeFarmSubset function with Farm*
{
	farmIDpresent(const std::vector<int> id_list) : ids(id_list) {} // constructor
	bool operator() (const Farm* f){ // overload operator function
		return (isWithin(f->Farm::get_id(),ids));
	}
	private:
		std::vector<int> ids; // member
};

#endif
