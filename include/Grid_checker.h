//© 2024 Colorado State University
#ifndef Grid_checker_h
#define Grid_checker_h

#include "Grid_cell.h"
#include "shared_functions.h"
#include "Local_spread.h"
#include "Status_manager.h"

extern int verboseLevel;

/// Makes comparisons between a focal farm and comparison farms in a cell
class Grid_checker
{
	private:
		const Parameters* p;
		int verbose; ///< Can be set to override global setting for console output
		std::vector<Grid_cell*> susceptible; ///< Local copy of cells, with vectors of susceptible farms within
//		std::vector<Farm*> exposed; ///< List of farms most recently exposed
		const std::unordered_map<int, Grid_cell*>* allCells; ///< Pointer to Grid_manager cells, referenced in infection evaluation among cells
        // variables for infection evaluation
        std::vector<std::string> speciesOnPrems; ///< List of species on all farms provided in premises file
        std::unordered_map<std::string,double> infExponents; ///< Species-specific infectiousness exponents, in same order as speciesOnAllFarms
		Status_manager* statusManagerPointer;
		Local_spread* kernel;
        int partial;
        std::vector<double> partialParams;
        std::tuple<double, double> latencyParams;

		void binomialEval(Farm* f1, Grid_cell* fc, Grid_cell* c2, int ccID, int t, int quarter_idx, std::vector<double> partialParams, std::vector<Farm*>& output, std::vector<double>& outputP); ///< Evaluates transmission from a focal farm to all susceptible farms in a cell via binomial method
//		void countdownEval(Farm*,Grid_cell*,Grid_cell*,int,std::vector<Farm*>&, int t, std::vector<double>partialParams); ///< Evaluates transmission from a focal farm to all susceptible farms in a cell via Keeling's "countdown" method
//		void pairwise(Farm*,Grid_cell*,Grid_cell*,int,std::vector<Farm*>&, int t , std::vector<double> partialParams); ///< Evaluates transmission from a focal farm to all susceptible farms in a cell pairwise

	public:
		///< Makes local copy of all Grid_cells, initially set as susceptible to check local spread against
		Grid_checker(const std::unordered_map<int, Grid_cell*>*,
			Status_manager*, const Parameters*);
		~Grid_checker();

		///< Function that handles actual comparisons between focal and susceptible premises.
		void stepThroughCells(
			std::vector<Farm*>&, // infectious
			std::vector<Farm*>&,//non-susceptible
            int t, int quarter_idx);
};

///> Determines if an object is present in a vector of objects.
/// Used to remove any Grid_cells that no longer contain any susceptible farms.
template<typename T>
struct isInList // used in grid checker for Grid_cell*
{
	isInList(const std::vector<T> in_list) : itemList(in_list) {} // constructor
	bool operator() (const T item){ // overload operator function
		return (isWithin(item,itemList));
	}
	private:
		std::vector<T> itemList; // member
};

#endif
