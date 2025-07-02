//© 2024 Colorado State University
#include "Grid_checker.h"

/// Makes shallow copy of Grid_cells to start as susceptible. Only the vector of pointers
/// to Farms is modified, not the Farms themselves (hence the shallow copy). Statuses are
/// only actually changed in Status_manager.
Grid_checker::Grid_checker(const std::unordered_map<int, Grid_cell*>* in_allCells,
	Status_manager* in_sm,
	const Parameters* p_in)
	:
	p(p_in),
	allCells(in_allCells),
    speciesOnPrems(p->species),
    infExponents(p->infExponents),
	statusManagerPointer(in_sm),
	kernel(p->kernel),
    partial(p->partial),
    partialParams(p->partialParams),
    latencyParams(p->latencyParams)

{
	verbose = verboseLevel;//verboseLevel;
	int fcount = 0;
	// initially copy all cells (containing all farms) into susceptible
	susceptible.reserve(allCells->size());
	for (auto& c:(*allCells)){
		susceptible.emplace_back(new Grid_cell(*(c.second))); // copy c and make new pointer to copy
		fcount += c.second->get_num_farms();
	}
	std::sort(susceptible.begin(),susceptible.end(),sortByID<Grid_cell*>);

if (verbose>1){std::cout<<"Grid checker constructed. "<<fcount<<" initially susceptible farms in "
	<<susceptible.size()<<" cells."<<std::endl;}
}

Grid_checker::~Grid_checker()
{
	for (auto& s:susceptible){delete s;}
}

/// Updates static list of cells with susceptible premises within. After transmission
/// evaluation, records sources of exposure in "sources" map from Status_manager and
/// records exposures in "exposed" vector, later accessed by Status_manager
/// \param[in] focalFarms	All currently infectious premises
/// \param[in] nonSus		All currently non-susceptible premises that cannot become exposed, including focalFarms
void Grid_checker::stepThroughCells(std::vector<Farm*>& focalFarms,
	std::vector<Farm*>& nonSus, int t, int quarter_idx)
{
//================================= update (vector of) susceptible Grid_cells*
	if (nonSus.size()>0){
if (verbose>1){std::cout<<"Removing "<<nonSus.size()<<" non-susceptible farms by ID."<<std::endl;}
		std::unordered_map<int, std::vector<int>> nonSusCellMap;
		for (auto& ns:nonSus){ // ns is Farm*
			nonSusCellMap[ns->Farm::get_cellID()].emplace_back(ns->Farm::get_id());
		}
		std::vector<Grid_cell*> empties;
		for (auto& s:susceptible){
			int sID = s->get_id(); // cell id
			if (nonSusCellMap.count(sID)==1){ // if this cell is represented in nonSus
if (verbose>1){std::cout<<"Removing matches in cell "<<sID<<std::endl;}
				s->removeFarmSubset(nonSusCellMap.at(sID)); // remove those farms from sus cell
				if (s->get_num_farms()==0){empties.emplace_back(s);}
			}
		}
		// remove any cells that have no (susceptible) farms left
		if (empties.size()>0){
			auto newSEnd = std::remove_if(susceptible.begin(),susceptible.end(),isInList<Grid_cell*>(empties));
if (verbose>1){
	int cCount = std::distance(newSEnd,susceptible.end());
	std::cout<<"Removing "<<cCount<<" cell(s) no longer susceptible."<<std::endl;
}
			susceptible.erase(newSEnd, susceptible.end());
		}
if (verbose>1){
	int scount = 0;
	for (auto& s:susceptible){scount += s->get_num_farms();}
	std::cout<<"Grid_checker:: Susceptible cells updated, now contain "<<scount<<" farms."<<std::endl;}

 	}

//================================= loop through pairs of inf farms and sus Grid_cells
	  // for each focal farm
	  for (auto& f1:focalFarms){
	  	int fcID = f1->Farm::get_cellID();
	  	Grid_cell* fc = allCells->at(fcID);
if (verbose>2){std::cout<<"Focal farm "<<f1->Farm::get_id()<<" in cell "<<fcID<<std::endl;}
		for (auto& c2:susceptible){
			int ccID = c2->Grid_cell::get_id();
			if (fc->kernelTo(ccID)>0){ // check if cell-cell tx possible
if (verbose>2){std::cout<<"Checking in-range comparison cell "<<ccID<<std::endl;}
			// Evaluation via gridding
			std::vector<Farm*> fToCellExp; // farms exposed by f1
			std::vector<double> trueProbs; // the true probability of transmission of each respective farm in fToCellExp
			binomialEval(f1, fc, c2, ccID, t, quarter_idx, partialParams, fToCellExp, trueProbs);
			// record sources of infection
			for (size_t exp_farm_idx=0; exp_farm_idx<fToCellExp.size(); ++exp_farm_idx){
                auto& exp1=fToCellExp.at(exp_farm_idx);
                double trueP=trueProbs.at(exp_farm_idx);
				statusManagerPointer -> Status_manager::add_premForEval(exp1, f1, 0, trueP);
			}

			//I don't think the following is used anymore, just leftovers since the pre-status_manager era /S
			// add to "grand total" exposed list from other farm-cell comparisons if not already present
//			for (auto& exp2:fToCellExp){
//				if (!isWithin<Farm*>(exp2,exposed)){
//					exposed.emplace_back(exp2);
//				}
//			}

			} // end if fc can infect cc
		} // end for loop through comparison cells
	  } // end for each focal farm

}

///	Calculates pmax of cell and N, draws h successes from binomial distribution
///	Randomly selects h farms, evaluates adjusted probabilities
/// \param[in]	f1	Infectious farm from which to evaluate transmission
///	\param[in]	fc	Focal cell containing infectious premises
///	\param[in]	c2	Comparison cell containing susceptible premises (can be same as fc)
///	\param[in]	ccID	ID of comparison cell
/// \param[in]  t  Timestep
/// \param[in]  partialParams  Parameters for the partial transmission function
///	\param[out] output  Vector of Farm*s exposed by this infectious farm
/// \param[out] outputP Vector of true probabilities of infection for each respective farm in the output vector.
void Grid_checker::binomialEval(Farm* f1, Grid_cell* fc, Grid_cell* c2, int ccID,
                                int t, int quarter_idx, std::vector<double> partialParams,
                                std::vector<Farm*>& output, std::vector<double>& outputP)
{


    double focalInfMax = f1->Farm::get_inf_max(); // the maximum infectiousness of a farm

    double focalInf; //current infectious of a farm
    //if the flag is 0 then use get_inf() to get the total infectiousness of farm
    if(p->infectionType == InfectionType::BTB)
    {
        int fid = f1->Farm::get_id();
        Prem_status* f1_pst = statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object
        focalInf = f1_pst->btb_get_prevalence_infectious();
    }
    else
    {
        if(partial==0){
            focalInf = f1->Farm::get_inf(quarter_idx);
        }else if(partial!=0){     //if the flag is 1 then use get_inf_partial_as_unvaccinated() to get the total infectiousness of farm as if it were unvaccinated
            int fid = f1->Farm::get_id();
            Prem_status* f1_pst = statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object
            focalInf = f1_pst->get_inf_partial_as_unvaccinated(t, quarter_idx, p, statusManagerPointer->get_normInf_map()); //the current infectiouness of a farm
        }else{
            std::cout<<"ERROR: In Grid_checker:: Partial tranisition flag does not exist. Exiting...";
            exit(EXIT_FAILURE);
        }
    }



	double kern = fc->Grid_cell::kernelTo(ccID); //if dangerousContacts_on, this includes DC prob
	double pmax; // Overestimated probability for any single premises
	if(p->infectionType == InfectionType::BTB)
    {
        pmax = focalInfMax * kern; //The kernel is already defined at a monthly timestep, so no need to scale up from daily timesteps.
    }
    else
    {
        pmax = oneMinusExp(-focalInfMax * kern * p->days_per_timestep);
    }

	double N = c2->Grid_cell::get_num_farms();
	std::vector<Farm*> fcexp; fcexp.reserve(N); // fcexp = "focal-comparison exposures"
	std::vector<double> fcexpP; fcexpP.reserve(N); // fcexpK = "focal-comparison exposures kernel value for (true) probabilities"

	// draw number of hypothetical farms exposed, from binomial
	int numExp = draw_binom(N, pmax);

	if (numExp == 0){ // no infected (or DC, if dangerousContacts_on) premises in this cell
	} else if (numExp > 0){
		// randomly choose numExp farms
		std::vector<Farm*> hypExposed; // hypothetically exposed
		std::vector<Farm*> compFarms = c2->get_farms();
		random_unique(compFarms,numExp,hypExposed);
		// evaluate each of the randomly selected farms
if(verbose>2){std::cout<<"Pmax: "<<pmax<<", "<<hypExposed.size()<<" hypothetical infections out of "
	<<compFarms.size()<<" farms in cell."<<std::endl;}
		std::vector<bool> dcEvaluations(p->NUM_DISEASE_STATUSES, false);
		for (auto& f2:hypExposed){
			// calc actual probabilities
			double f1x = f1 -> Farm::get_x();
			double f1y = f1 -> Farm::get_y();
			double f2x = f2 -> Farm::get_x();
			double f2y = f2 -> Farm::get_y();
			double xdiff = (f1x - f2x);
			double ydiff = (f1y - f2y);
			double distBWfarmssq = xdiff*xdiff + ydiff*ydiff;

			// calculate probability between these specific farms
			double ptrue;
			if(p->infectionType == InfectionType::BTB)
            {
                double kernelBWfarms = kernel->atDistSq(distBWfarmssq, f1->get_parent_county()->get_wildlife_density(), quarter_idx); // kernelsq calculates kernel based on distance squared
                ptrue = focalInf * kernelBWfarms; //The kernel is already defined at a monthly timestep, so no need to scale up from daily timesteps.
                if(f1 == f2)
                {
                    ptrue = 0.0; //For btb, premises *can* be exposed multiple times through local spread, but a prem can't expose itself.
                }
            }
			else
            {
                double kernelBWfarms = kernel->atDistSq(distBWfarmssq); // kernelsq calculates kernel based on distance squared
                double compSus = f2->Farm::get_sus(quarter_idx); // susceptible farm in comparison cell
                ptrue = oneMinusExp(-focalInf * compSus * kernelBWfarms * p->days_per_timestep); // prob tx between this farm pair
                if(verbose>2){std::cout<<"Inf: "<<focalInf<<", sus: "<<compSus<<", kernel: "<<kernelBWfarms<<", Ptrue "<<ptrue<<std::endl;}
            }

			double random = uniform_rand();
			if (random <= ptrue/pmax){ // actual infection
if (verbose>1){std::cout << "Infection @ distance: "<< std::sqrt(distBWfarmssq)/1000 << " km, prob "<<ptrue<<std::endl;}
				fcexp.push_back(f2);
				fcexpP.push_back(ptrue);
			}
			std::string f1Status = statusManagerPointer -> Status_manager::getAny_fileStatus(f1);
			if (p->dangerousContacts_on==1 && f1Status.compare("reported")!=0){
			// only bother calculating DCs if f1 is not yet reported, any DCs post-reporting are never used
				std::fill(dcEvaluations.begin(), dcEvaluations.end(), false);
				bool possibleDC = 0;
				for (auto& r:(p->dcRiskScale)){
					double pDC = ptrue*r.second; // prob of being DC when status is r.first
					double random = uniform_rand();
					int status_idx = static_cast<int>(r.first); //DiseaseStatus to int.
					if (random <= pDC/pmax){
						dcEvaluations[status_idx] = true;
						possibleDC = 1;
if (verbose>1){std::cout<<"Dangerous contact identified"<<std::endl;}
					} else {
						dcEvaluations[status_idx] = false;
					}
				}
				if (possibleDC){
					// store DC evaluations with f1
					if(verbose>1){std::cout<<"GC before call to status manager add_potentialDC"<<std::endl;}
					statusManagerPointer -> Status_manager::add_potentialDC(f1, f2, dcEvaluations);
					if(verbose>1){std::cout<<"GC after call to status manager add_potentialDC"<<std::endl;}
				}
			}

		 } // end "for each hypothetically exposed farm"
		} // end "if any hypothetically exposed farms"
	fcexp.swap(output);
	fcexpP.swap(outputP);
}

//Unused local speread transmission function (Keelings method)
/////	Calculates and evaluates probability of cell entry, then steps through each premises
/////	in cell and evaluates adjusted probabilities (Keeling's method)
///// \param[in]	f1	Infectious farm from which to evaluate transmission
/////	\param[in]	fc	Focal cell containing infectious premises
/////	\param[in]	c2	Comparison cell containing susceptible premises (can be same as fc)
/////	\param[in]	ccID	ID of comparison cell
/////	\param[out]	output	Vector of Farm*s exposed by this infectious premises
//void Grid_checker::countdownEval(Farm* f1, Grid_cell* fc, Grid_cell* c2, int ccID,
//	std::vector<Farm*>& output, int t,  std::vector<double> partialParams)
//{
//    double focalInfMax = f1->Farm::get_inf_max(); // the maximum infectiousness of a farm
//
//    double focalInf; //current infectious of a farm
//    //if the flag is 0 then use get_inf() to get the total infectiousness of farm
//    if(partial==0){
//        focalInf = f1->Farm::get_inf();
//    }else if(partial!=0){     //if the flag is 1 then use get_inf_partial_as_unvaccinated() to get the total infectiousness of farm
//        int fid = f1->Farm::get_id();
//        Prem_status* f1_pst = statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object
//        focalInf = f1_pst->get_inf_partial_as_unvaccinated(t, p, statusManagerPointer->get_normInf_map()); //the current infectiouness of a farm
//    }else{
//        std::cout<<"ERROR: In Grid_checker:: Partial tranisition flag does not exist. Exiting...";
//        exit(EXIT_FAILURE);
//    }
//	double kern = fc->Grid_cell::kernelTo(ccID);
//	double pmax = oneMinusExp(-focalInfMax * kern * p->days_per_timestep); // Overestimated probability for any single premises, "prob6" in MT's Fortran code:
//	double N = c2->Grid_cell::get_num_farms();
//	double pcell = oneMinusExp(-focalInf * kern * N * p->days_per_timestep); // Probability of cell entry
//
//	std::vector<Farm*> fcexp; fcexp.reserve(N); // fcexp = "focal-comparison exposures"
//
//	double s = 1; // on/off switch, 1 = on (single hypothetical infection hasn't happened yet)
//	double random1 = uniform_rand();
//// Grid checkpoint A
//	if (random1 <= pcell){ // if farm to cell succeeds
// 		int f2count = 1; // how many farms in comparison cell have been checked
// 		std::vector<Farm*> compFarms = c2->get_farms();
//		for (auto& f2:compFarms){
//			double oneMinusExpA = oneMinusExp(-focalInf * kern * (N+1-f2count) * p->days_per_timestep); // 1 - exp(A)
//			double pcellAdj = 1 - s + s*oneMinusExpA; // equivalent to 1 - s*exp(A)
//			double random2 = uniform_rand(); // "prob4" in MT's Fortran code
//// Grid checkpoint B
//			if (random2 <= pmax/pcellAdj){
//			// if (one max susceptible)/(entrance prob accounting for # of farms checked) succeeds
//			s = 0; // remainingFarmProb recalculates to 1 for remainder of loop
//			// get actual distances between farms
//			double f1x = f1 -> Farm::get_x();
//			double f1y = f1 -> Farm::get_y();
//			double f2x = f2 -> Farm::get_x();
//			double f2y = f2 -> Farm::get_y();
//			double xdiff = (f1x - f2x);
//			double ydiff = (f1y - f2y);
//			double distBWfarmssq = xdiff*xdiff + ydiff*ydiff;
//			double kernelBWfarms = kernel->atDistSq(distBWfarmssq); // kernelsq calculates kernel based on distance squared
//			double compSus = f2->Farm::get_sus(); // susceptible farm in comparison cell (farmInf already defined from focal cell)
//
//			// calculate probability between these specific farms
//			// "prob3" in MT's Fortran code
//			double ptrue = 1-exp(-focalInf * compSus * kernelBWfarms);
//// Grid checkpoint C
//			if (random2 <= ptrue/pcellAdj){
//				// infect
//				if(verbose>1){
//					std::cout << "Infection @ distance: ";
//					std::cout << std::sqrt(distBWfarmssq)/1000 << ", prob "<<ptrue<<std::endl;
//				}
//				fcexp.emplace_back(f2);
//			}
//		 } // end "if farm hypothetically exposed"
//		 f2count++;
//		} // end "for each comparison farm"
//	} // end "if >1 hypothetical infection"
//	fcexp.swap(output);
//}


//Unused unoptimzed naive pairwise local spread evaluation.
///// Calculates filtered pairwise transmission: only makes pairwise calculations if random
///// number passes pmax filter
///// \param[in]	f1	Infectious farm from which to evaluate transmission
/////	\param[in]	fc	Focal cell containing infectious premises
/////	\param[in]	c2	Comparison cell containing susceptible premises (can be same as fc)
/////	\param[in]	ccID	ID of comparison cell
/////	\param[out]	output	Vector of Farm*s exposed by this infectious premises
//void Grid_checker::pairwise(Farm* f1, Grid_cell* fc, Grid_cell* c2, int ccID,
//	std::vector<Farm*>& output, int t,  std::vector<double> partialParams)
//{
//    double focalInfMax = f1->Farm::get_inf_max(); // the maximum infectiousness of a farm
//
//    double focalInf; //current infectious of a farm
//    //if the flag is 0 then use get_inf() to get the total infectiousness of farm
//    if(partial==0){
//        focalInf = f1->Farm::get_inf();
//    }else if(partial!=0){     //if the flag is 1 then use get_inf_partial_as_unvaccinated() to get the total infectiousness of farm
//        int fid = f1->Farm::get_id();
//        Prem_status* f1_pst = statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object
//        focalInf = f1_pst->get_inf_partial_as_unvaccinated(t, p, statusManagerPointer->get_normInf_map()); //the current infectiouness of a farm
//    }else{
//        std::cout<<"ERROR: In Grid_checker:: Partial tranisition flag does not exist. Exiting...";
//        exit(EXIT_FAILURE);
//    }
//
//	double kern = fc->Grid_cell::kernelTo(ccID);
//	double pmax = oneMinusExp(-focalInfMax * kern * p->days_per_timestep); // Overestimate of p for all farms in this cell
//	std::vector<Farm*> cFarms = c2->Grid_cell::get_farms();
//
//	std::vector<Farm*> fcexp; fcexp.reserve(cFarms.size()); // For output
//
//	for (auto& cf:cFarms){
//		double random = uniform_rand();
//		if (random <= pmax){
//			double f1x = f1 -> Farm::get_x();
//			double f1y = f1 -> Farm::get_y();
//			double f2x = cf -> Farm::get_x();
//			double f2y = cf -> Farm::get_y();
//			double xdiff = (f1x - f2x);
//			double ydiff = (f1y - f2y);
//			double distBWfarmssq = xdiff*xdiff + ydiff*ydiff;
//			double kernelBWfarms = kernel->atDistSq(distBWfarmssq); // kernelsq calculates kernel based on distance squared
//			double compSus = cf->Farm::get_sus(); // susceptible farm in comparison cell
//
//			// calculate probability between these specific farms
//			double ptrue = oneMinusExp(-focalInf * compSus * kernelBWfarms * p->days_per_timestep); // prob tx between this farm pair
//			if (random <= ptrue){ // actual infection
//if(verbose>1){std::cout << "Infection @ distance: "<< std::sqrt(distBWfarmssq)/1000 << " km, prob "<<ptrue<<std::endl;}
//					fcexp.emplace_back(cf);
//			}
//		}
//	}
//	fcexp.swap(output);
//}
