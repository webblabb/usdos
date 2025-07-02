//© 2024 Colorado State University
#include "Status_manager.h"
#include "shared_functions.h"

#include <algorithm>

/// Establishes sequences of statuses for disease, file status, and control statuses.
/// Sets seed farm(s) as exposed, using index reporting time for reporting.
Status_manager::Status_manager(std::vector<Farm*>& focalFarms, const Parameters* parameters,
	Grid_manager* grid, Control_manager* control, Diagnostic_manager* diagnostic) :
		seededFarms(focalFarms), // saved for output
		parameters(parameters),
		controlManager(control),
		diagnosticManager(diagnostic),
		gridManager(grid),
    allPrems(grid->get_allFarms()),
    allCounties(grid->get_allCountiesMap()),
    allStates(grid->get_allStates()),
    allControlTypes(control->get_controlTypes()),
    controlResources(control->get_controlResources()),
    controlReleaseSchedule(control->get_resourceBoostSchedule()),
    allDiagnosticTypes(diagnostic->get_diagnosticTypes()),
    diagnosticResources(diagnostic->get_diagnosticResources()),
    recentNotSus(0),
    pastEndTime(std::make_tuple(parameters->timesteps+100, 0)),
    nPrems(allPrems->size()),
    species(parameters->species) // store species for formatting later
{
	verbose = verboseLevel;
	btb_slaughter_ofn = parameters->batchDateTime + "_slaughter.txt";


	// Specify duration of each disease status and what follows
	if(parameters->infectionType == InfectionType::FMD)
    {
        statusShift exp {parameters->latencyParams, "inf"};
        statusSequences["exp"] = exp;

        statusShift inf {parameters->infectiousParams, "imm"};
        statusSequences["inf"] = inf;

        statusShift imm {pastEndTime, "NA"};
        statusSequences["imm"] = imm;
    }
    else if(parameters->infectionType == InfectionType::BTB)
    {
        statusShift btb {pastEndTime, "rec"}; //Affected by btb, doesn't enxpire after a set amount of time, but when animals die off or are sent from the prem.
        statusSequences["btb"] = btb;

        statusShift rec {pastEndTime, "btb"}; //BTB ONLY. Recovered from btb, last infected animals were removed and the prem is back to being treated as susceptible.
        statusSequences["rec"] = rec;
    }

	// set seed farms as exposed
	for (auto& f:focalFarms)
	{
        int fid = f->Farm::get_id();

        if(parameters->infectionType == InfectionType::FMD)
        {
            set_status(f, 1, "exp", parameters->latencyParams); // "exp" = disease status exposed
        }
        else if(parameters->infectionType == InfectionType::BTB)
        {
            set_status(f, 1, "btb", pastEndTime);
            Prem_status* seeded_pst = changedStatus.at(fid);
            btb_n_freshly_infected_animals += seeded_pst->btb_seed();
            //Guarantees that the seeded prem is recorded as infected, even if it recovers during the first time step.
            //Won't save duplicate records since btb_all_infected is a set.
            btb_all_infected.insert(seeded_pst);
        }
        else
        {
            std::cout << "Unknown type of infection simulated." << std::endl;
            exit(EXIT_FAILURE);
        }

        //set exposure time for source farms
        //record exposedBy by itself
        //route is local spread
        //not blocked
        changedStatus.at(fid)->Prem_status::add_exposureSource(1, f, 0, "none");
	}

	//the status shift sequence is only for things that are controlled only by a time delay
	if (parameters->control_on == true || parameters->diagnostic_on == true){
			if (parameters-> testSuspects_on == false){
				statusShift exposed {pastEndTime, "reported"};
			// End times for fileStatus 'exposed' differs depending on dangerousContact status, handled in expose() and set_status(). pastEndTime put in as placeholder.
			statusSequences["exposed"] = exposed;
			statusShift reported {pastEndTime, "NA"}; // 'reported' status is permanent
			statusSequences["reported"] = reported;
		 }else if (parameters-> testSuspects_on == true){
			statusShift exposed {pastEndTime, "suspected"};
			// End times for fileStatus 'exposed' differs depending on dangerousContact status, handled in expose() and set_status(). pastEndTime put in as placeholder.
			statusSequences["exposed"] = exposed;
			statusShift suspected {pastEndTime, "NA"}; // 'reported' status is permanent
			statusSequences["suspected"] = suspected;
			}
		}
		if (parameters->diagnostic_on == true){
			for (auto& dt:(parameters->diagnosticTypes)){
				std::string d_type = dt;
				//statusShift exposed {pastEndTime, "suspected"};
				// End times for fileStatus 'exposed' differs depending on dangerousContact status, handled in expose() and set_status(). pastEndTime put in as placeholder.

				statusShift started {parameters->testStartToCompleteLag.at(d_type), "complete."+d_type};
				statusSequences["started."+d_type] = started;

				statusShift complete {pastEndTime, "reported"};
				statusSequences["complete."+d_type] = complete;

				statusShift reported {pastEndTime, "NA"};
				statusSequences["reported"] = reported;
			}
		}

			// Specify duration of each control status and what follows
		if (parameters->control_on == true){
			for (auto& ct:(parameters->controlTypes)){
				std::string c_type = ct;

				statusShift implemented {parameters->implementToEffectiveLag.at(c_type), "effective."+c_type};
				statusSequences["implemented."+c_type] = implemented;

				statusShift effective {parameters->effectiveToInactiveLag.at(c_type), "inactive."+c_type};
				statusSequences["effective."+c_type] = effective;

				statusShift inactive {pastEndTime, "NA"};
				statusSequences["inactive."+c_type] = inactive;

			}
		}

        if (parameters->testSuspects_on == false){
            // set seed farms as exposed (fileStatus), for reporting
            std::tuple<double, double> iReport = parameters->indexReportLag;
            for (auto& f:focalFarms){
                int indexReportTime = normDelay(iReport); // current time + index report delay
                std::tuple<double, double> indexReportTuple = std::make_tuple(indexReportTime, 0);
                set_status(f, 1, "exposed", indexReportTuple); // file status
            }
        }else if (parameters->testSuspects_on == true){
                // set seed farms as exposed (fileStatus), for suspecting
                std::tuple<double, double> iSuspect = parameters->indexSuspectLag;
                for (auto& f:focalFarms){
                    int indexSuspectTime = normDelay(iSuspect); // current time + index report delay
                    std::tuple<double, double> indexSuspectTuple = std::make_tuple(indexSuspectTime, 0);
                    set_status(f, 1, "exposed", indexSuspectTuple); // file status
            }
        }

//if (verbose>1){
//	std::cout<<focalFarms.size()<<" exposed premises initiated. End times for exposed premises: "<<std::endl;
//	for (auto& e : diseaseStatuses.at("exp").units){
//		std::cout << e->get_end("exp") << ", ";
//	}
//	std::cout<<std::endl;
//}

if(verbose > 0){
	std::cout << "Status manager initiated."<<std::endl;
	}
}

Status_manager::~Status_manager()
{
    for (auto& f:changedStatus){delete f.second;}
	for (auto& c:changedCoStatus){delete c.second;}
	for (auto& c:changedStateStatus){delete c.second;}
	if(btb_slaughter_ofs.is_open())
    {
        btb_slaughter_ofs.close();
    }
}

/// Checks if a Prem_status exists for a premises and creates one if not. Returns
/// identifier for premises that is used as index value in changedStatus.
int Status_manager::verify_premStatus(Farm* f){
	int fid = f->Farm::get_id();
	// if a prem_status object has not yet been made for this farm, make one
	if (changedStatus.count(fid)==0){
		changedStatus[fid] = new Prem_status(f);
	}
	return fid;
}

/// Records confirmed exposure sources for a Prem_status without having to set disease
/// status (i.e. allows tracking of exposures avoided by vaccination, where disease status
/// is still susceptible). Used by filter_shipments() for shipments that are
/// banned, and eval_exposure() for realized premises-level control
void Status_manager::add_premSource(int t, Farm* toBeExposed, Farm* exposedBy,
	int route, std::string prevented)
{
	int fid = verify_premStatus(toBeExposed);
	changedStatus.at(fid)->Prem_status::add_exposureSource(t, exposedBy, route, prevented);
	// Also add info to "sources" for printing
	sources.emplace_back(std::make_tuple(toBeExposed, exposedBy, route, prevented));
}

/// Changes status of a farm, and determines time when next status begins.
/// Adds farm to appropriate status list with start and end times. If this is the first
/// farm listed for a particular status, sets the placeholder for that status at the
/// beginning of the status list.
/// \param[in]	f	Farm for which status should be changed
///	\param[in]	startTime	Time at which this status begins
/// \param[in]	status status to which this farm will be set
///	\param[in]	statusDuration	Tuple of mean, variance in days until next status begins. For permanent statuses (i.e. reported), use mean=pastEndTime, var=0
/// \param[in] 	controlEffect Optional argument needed when setting control to effective or inactive

void Status_manager::set_status(Farm* f, int startTime, std::string status,
	std::tuple<double,double> statusDuration, std::tuple<double, double> controlEffect)
{

	int fid = verify_premStatus(f);
	Prem_status* pst = changedStatus.at(fid);

    if(status == "started")
    {
        std::cout << "sdd" << std::endl;
    }
	// if lag time = 0, use next status & duration instead
	while (std::get<0>(statusDuration) == 0){
		if (statusSequences.at(status).next.compare("NA")!=0){ // if there is a next status
			status = statusSequences.at(status).next; // overwrite status as the next status
			statusDuration = statusSequences.at(status).duration; // get duration of that status
		} else {
			// leave status as is and set duration to pastEndTime
			statusDuration = pastEndTime;
		}
	}

if(verbose>1){std::cout<<"SM::set_status: Setting "<<fid<<" status to ";}

	// if status is a disease status:
	if (status.compare("exp")==0 ||
        status.compare("inf")==0 ||
        status.compare("imm")==0)
    {
		if(verbose>1){std::cout<<status<<" (disease status)"<<std::endl;}

		changedStatus.at(fid)->Prem_status::set_diseaseStatus(status);
		// add to appropriate statusList
		if (diseaseStatuses.count(status)==0){
			diseaseStatuses[status].lo = 0; // set placeholder at beginning
		}
        diseaseStatuses.at(status).units.emplace_back(changedStatus.at(fid));

		if (status.compare("exp")==0){ // if farm is becoming exposed
			notSus.emplace_back(allPrems->at(fid)); // add to list of non-susceptibles
		}

	} else if (status.compare("btb")==0 ||
               status.compare("rec")==0)
    {
	    pst->Prem_status::set_diseaseStatus(status);
		// add to appropriate statusList
		if (diseaseStatuses.count(status)==0){
			diseaseStatuses[status].lo = 0; // set placeholder at beginning
		}
		//Premises can be exposed multiple times with bTB, but only one prem status
		//pointer should be added here even if the prem is exposed multiple times.
		//This doesn't happen for FMD since an exposed/inf prem can't be infected again,
		//so check if any of the pointers in diseaseStatuses.at(status) is equal to
		//this premises' prem_status pointer.
		if(std::none_of(diseaseStatuses.at(status).units.begin(),
                        diseaseStatuses.at(status).units.end(),
                        [pst](Prem_status* pst_comp) { return pst==pst_comp; }))
        {
            diseaseStatuses.at(status).units.emplace_back(pst);
        }

        //When simulating btb, we don't want premises to be removed from the premises that potentially
        //can be exposed just because they get exposed once - they can be exposed multiple times.
//		if (status.compare("btb") == 0){ // if farm is becoming exposed
//			notSus.emplace_back(allPrems->at(fid)); // add to list of non-susceptibles
//		}
    // if status is a file status (applies to all types of control):
	// other fileStatus: notDangerousContact
	}  if (status.compare("dangerousContact")==0 ||
           status.compare("exposed")==0 ||
           status.compare("suspected")==0 ||
           status.compare("reported")==0)
    {
		if(verbose>1)
			{std::cout<<status<<" (file status)"<<std::endl;}

		changedStatus.at(fid)->Prem_status::set_fileStatus(status);
		// add to appropriate file-status list
		if (fileStatuses.count(status)==0){
			fileStatuses[status].lo = 0; // set placeholder at beginning
		}
		fileStatuses.at(status).units.emplace_back(changedStatus.at(fid));

		if (status.compare("suspected")==0){
			newPremSuspects.emplace_back(changedStatus.at(fid));
		}
		if (status.compare("reported")==0){
			newPremReports.emplace_back(changedStatus.at(fid));
			if (parameters->dangerousContacts_on==1){
				// determine which of potential dangerous contacts will be dangerous contacts
				std::unordered_map<Farm*, std::vector<bool>> pDCs =
					*(changedStatus.at(fid)->Prem_status::get_potentialDCs());
if(verbose>1 && pDCs.size()>0){
	std::cout<<"SM::Reported farm has "<<pDCs.size()<<" potential DCs"<<std::endl;}
				for (auto& dc:pDCs){ // dc.first is prem to evaluated as DC
if (verbose>1)std::cout<<"SM::set_status inside loop over DCs"<<std::endl;
if (verbose>1){std::cout<<dc.first->Farm::get_id()<<" is DC/reported? "<<getAny_fileStatus(dc.first)<<std::endl;}
					if (getAny_fileStatus(dc.first).compare("dangerousContact")!=0 ||
						//getAny_fileStatus(dc.first).compare("suspected")!=0 ||
						getAny_fileStatus(dc.first).compare("reported")!=0){ // if not already suspected, reported or already a DC

						// get current disease status
						std::string dcDiseaseStatusStr = getAny_diseaseStatus(dc.first);
						DiseaseStatus dcDiseaseStatus = parameters->StringToDiseaseStatus.at(dcDiseaseStatusStr);
						int status_idx = static_cast<int>(dcDiseaseStatus);

						// determine if dc should be dangerousContact
						if (dc.second.at(status_idx)){ // dc.second is a vector of bools
							bool isDC = dc.second.at(status_idx);
	if(verbose>1){std::cout<<"SM::set_status: potential DC "<<dc.first->Farm::get_id()<<" disease status = "<<dcDiseaseStatusStr<<", isDC = "<<isDC<<std::endl;}
							bool sourceIsF = 1; // assume source is reported farm unless changed...
							if (dcDiseaseStatus != DiseaseStatus::SUS and dcDiseaseStatus != DiseaseStatus::REC){ // if not susceptible, was it exposed by a different farm?
								sourceIsF = changedStatus.at(dc.first->Farm::get_id())->Prem_status::is_exposureSource(f);
							}
	if(verbose>1){std::cout<<"SM::set_status: potential DC "<<dc.first->Farm::get_id()<<"'s source is reported farm? "<<sourceIsF<<std::endl;}
							if (isDC==1 && sourceIsF==1){
								changedStatus.at(fid) -> Prem_status::add_DC(dc.first);
								set_status(dc.first, startTime, "dangerousContact", parameters->dcReportLag);
							}

							} // end "if status has a corresponding DC outcome"
					} // end "if not already a DC or reported"

				} // end "for each potential dc"
				dcsPerIP.emplace_back((changedStatus.at(fid) -> Prem_status::get_dangerousContacts()).size());
			} // end "if dangerousContacts_on"
		} // end "if being reported now"
	}

    // if status is a diagnostic status (specific to each premises-level diagnostic type):
    if (// first n letters of status are:
        status.compare(0,7,"started")==0 ||
        status.compare(0,8,"complete")==0)
    {

        if(verbose>1)
        {
            std::cout<<status<<" (diagnostic status)"<<std::endl;
        }
        changedStatus.at(fid)->Prem_status::set_diagnosticStatus(status);


        // add to appropriate diagnostic-status list
        if (diagnosticStatuses.count(status)==0)
        {
            diagnosticStatuses[status].lo = 0; // set placeholder at beginning
        }
        diagnosticStatuses.at(status).units.emplace_back(changedStatus.at(fid));


        // get the diagnostic type at the end of the status string
        std::string d_type = status;
        d_type.replace(0, d_type.find(".")+1, "");
        //std::tuple<double, double> sensitivity = allDiagnosticTypes->at(d_type)->sensitivity;


        if(status.compare(0,8,"complete")==0)
        {

            int fid = f->Farm::get_id();
            //Ugly fix here. If Lag between testing and test complete is zero, as it likely will be for btb since the timestep is 1 month,
            //the actual started status will be skipped due to how this function is constructed (skips past all statuses with duration 0).
            //Since the code that actually performs the test happens after the "started" status of this prem has been skipped
            //it will never be executed if lag = 0, and no farm will be tested.
            //To get around that I've made this fix that if simulating btb and lag is 0 will do the test here instead, once
            //the test is considered completed - meaning during the same timestep that it got started. I think this works, assuming
            //that the test is never also performed outside of here.
            if(parameters->infectionType == InfectionType::BTB and
               std::get<0>(parameters->testStartToCompleteLag.at(d_type)) == 0)
            {
                bool testPositive = eval_probHerdPositive(f, -1, -1, d_type); //Second parameter is current timestep and third is current quarter, but they are only used inside eval_probHerdPositive if simulating FMD so can be set to -1 for bTB.
                pst->Prem_status::add_testResult(d_type, 0);

                if(testPositive == 0){
                    if (verbose>1){std::cout<<"Diagnostic test: negative"<<std::endl;}
                    pst->Prem_status::add_testResult(d_type, 0);
                }else if(testPositive == 1){
                    if (verbose>1){std::cout<<"Diagnostic test: positive"<<std::endl;}
                    pst->Prem_status::add_testResult(d_type, 1);
                }
            }

            //Here the function is resumed as usual, so if simulating FMD, or BTB with lag > 0, things should work as before.
            if (getAny_fileStatus(f).compare("reported") != 0)
            {
                int diagRes =  changedStatus.at(fid)->Prem_status::get_testResult(d_type);
                if (diagRes == 1)
                {
                    set_status(f,startTime,"reported", pastEndTime);
                }
                else if(diagRes != 1 and parameters->infectionType == InfectionType::BTB)
                {
                    set_status(f,startTime,"exposed", pastEndTime);
                }
            }
        }
    }

	// if status is a control status (specific to each premises-level control type):
	if (// first n letters of status are:
			status.compare(0,11,"implemented")==0 ||
			status.compare(0,9,"effective")==0 ||
			status.compare(0,8,"inactive")==0){

		if(verbose>1){std::cout<<status<<" (control status)"<<std::endl;}

		// add to appropriate file-status list
		if (controlStatuses.count(status)==0){
			controlStatuses[status].lo = 0; // set placeholder at beginning
		}
		controlStatuses.at(status).units.emplace_back(changedStatus.at(fid));

		// get the control type at the end of the status string
		std::string c_type = status;
		c_type.replace(0, c_type.find(".")+1, "");
		std::tuple<double, double> effectiveness = allControlTypes->at(c_type)->effectiveness;

		if (status.compare(0,9,"effective")==0) {
            changedStatus.at(fid)->Prem_status::add_controlStatus(c_type); // used in detailed output when tx/exp prevented
            if(c_type == "vax"){
				newPremVaxs.emplace_back(changedStatus.at(fid));
                changedStatus.at(fid)->vaccinate(std::get<0>(effectiveness)); //For vaccination "effectiveness" is the vaccine efficacy on the animal level and is stored in element 0 of the effectiveness tuple.
                changedStatus.at(fid)->Prem_status::add_probPreventExposure(0.0); //For vaccination the prevention is not applied as for prem-level control strategies.
                changedStatus.at(fid)->Prem_status::add_probPreventTransmission(0.0);
            } else{
                changedStatus.at(fid)->Prem_status::add_probPreventExposure(std::get<0>(effectiveness));
                changedStatus.at(fid)->Prem_status::add_probPreventTransmission(std::get<1>(effectiveness));
            }
		} else if (status.compare(0,8,"inactive")==0) {
		    if(c_type == "vax"){
		        changedStatus.at(fid)->unvaccinate();
				changedStatus.at(fid)->Prem_status::rem_probPreventExposure(0.0);
				changedStatus.at(fid)->Prem_status::rem_probPreventTransmission(0.0);
		    } else {
				changedStatus.at(fid)->Prem_status::rem_probPreventExposure(std::get<0>(effectiveness));
				changedStatus.at(fid)->Prem_status::rem_probPreventTransmission(std::get<1>(effectiveness));
		    }
            //If btb is being simulated, the prem returns to being susceptible once the control "wears off".
            if(parameters->infectionType == InfectionType::BTB){
                    pst->btb_remove_all_infected();
                    set_status(pst, startTime, "btb", std::make_tuple<double, double>(startTime+1, 0));
            }
		}
	}


	// set start and end times
	changedStatus.at(fid)->Prem_status::set_start(status,startTime);
	int endTime = startTime+normDelay(statusDuration);
	changedStatus.at(fid)->Prem_status::set_end(status,endTime);
}

void Status_manager::set_regionStatus(std::string id, std::string regionType, int startTime, std::string status,
	std::tuple<double,double> statusDuration, std::tuple<double,double> controlEffect)
{
	// if a Region_status object has not yet been made, make one
	Region_status* regionStatusPointer;
	if (regionType.compare("county")==0){
		if (changedCoStatus.count(id)==0){
			if (allCounties->count(id)<1){
				std::cout<<"ERROR in Status_manager::set_regionStatus: Attempting to set status for a county that does not exist in county map. Exiting...";
				exit(EXIT_FAILURE);
			}
			changedCoStatus[id] = new Region_status(allCounties->at(id));
		}
		regionStatusPointer = changedCoStatus.at(id);
	} else if (regionType.compare("state")==0){ // region is state
		if (changedStateStatus.count(id)==0){
			if (allStates->count(id)<1){
				std::cout<<"ERROR in Status_manager::set_regionStatus: Attempting to set status for a state that does not exist in state map. Exiting...";
				exit(EXIT_FAILURE);
			}
 			changedStateStatus[id] = new Region_status(allStates->at(id));
	 	}
	 	regionStatusPointer = changedStateStatus.at(id);
	}

	// if lag time = 0, use next status & duration instead
	while (std::get<0>(statusDuration) == 0){
		if (statusSequences.at(status).next.compare("NA")!=0){ // if there is a next status
			status = statusSequences.at(status).next; // overwrite status as the next status
			statusDuration = statusSequences.at(status).duration; // get duration of that status
		} else {
			// leave status as is and set duration to pastEndTime
			statusDuration = pastEndTime;
		}
	}

	if(verbose>1){std::cout<<"SM::set_regionStatus: Setting "<<regionType<<" "<<id<<" status to "<<status
	<<" for "<<std::get<0>(statusDuration)<<" days"<<std::endl;}

	// Disease statuses are not tracked at the region level, only file and control statuses

	if (// first n letters of status are:
		status.compare(0,11,"implemented")==0 ||
		status.compare(0,9,"effective")==0 ||
		status.compare(0,8,"inactive")==0){

		regionStatusPointer->Region_status::add_controlStatus(status); // not used for anything, but stored just in case
		// add to appropriate statusList
		if (regionControlStatuses.count(status)==0){ // if this is the first region for this status
			regionControlStatuses[status].lo = 0; // set placeholder at beginning
		}
		regionControlStatuses.at(status).units.emplace_back(regionStatusPointer);
		// get the control type at the end of the status string
		std::string c_type = status;
		c_type.replace(0, c_type.find(".")+1, "");
		std::tuple<double, double> effectiveness = allControlTypes->at(c_type)->effectiveness;

		if (status.compare(0,9,"effective")==0){
				regionStatusPointer->Region_status::add_probPreventExposure(std::get<0>(effectiveness));
				regionStatusPointer->Region_status::add_probPreventTransmission(std::get<1>(effectiveness));
		} else if (status.compare(0,8,"inactive")==0){
				regionStatusPointer->Region_status::rem_probPreventExposure(std::get<0>(effectiveness));
				regionStatusPointer->Region_status::rem_probPreventTransmission(std::get<1>(effectiveness));
		}
	}

	// set start and end times
	regionStatusPointer->Region_status::set_start(status,startTime);
	int endTime = startTime+normDelay(statusDuration);
	regionStatusPointer->Region_status::set_end(status,endTime);

}

/// "Removes" expired farms with a given status and if applicable, sets next status.
/// \param[in] t Current timestep, used to determine if statuses have expired and set next status
/// \param[in] inStatus statusList containing vector of farms and integer placeholder
/// \param[in] status Status that applies to inStatusList
void Status_manager::update(int t, std::string status, statusList<Prem_status*>& inStatusList)
{

    // check from iterator "lo" to end, switch expired farms to before "lo"
    if (inStatusList.lo < inStatusList.units.size())  // not all farms have expired
    {
        std::vector<Prem_status*>::iterator lo_it = inStatusList.units.begin();
        std::advance(lo_it, inStatusList.lo); // move iterator to lo

        for (auto it = lo_it; it < inStatusList.units.end(); it++)  // check each farm* between lo and end
        {
            if ( (*it)->Prem_status::get_end(status) == t )  // if validity of this status expires for this farm today
            {
                // set to next status
                if (status.compare("exposed")==0 && parameters->testSuspects_on == false)  // specific to file status, not disease
                {
                    // report and add to waitlist for each control type
                    //If we want to add a coin test to determine whether an exposed farm will be reported, then add here
                    if(get_totalPremsWithFileStatus("reported")<1)
                    {
                        firstRepTime = t; //if the number of reporteds is less than 1 then record the first reported time.

                        if(verbose>1)
                        {
                            std::cout << "SM:: First report time recorded as "<< firstRepTime << std::endl;
                        }
                    }
                    report_countyAndState(*it, t); // also report at county and state level
                }
                else if (status.compare(0,7,"started")==0 && parameters->testSuspects_on == true)    // specific to file status, not disease
                {
                    // report and add to waitlist for each control type
                    //If we want to add a coin test to determine whether an exposed farm will be reported, then add here
                    if(get_totalPremsWithFileStatus("reported")<1)
                    {
                        firstRepTime = t; //if the number of reporteds is less than 1 then record the first reported time.
                        if(verbose>1)
                        {
                            std::cout << "SM:: First report time recorded as "<< firstRepTime << std::endl;
                        }
                    }
                    report_countyAndState(*it, t); // also report at county and state level
                }

                if (statusSequences.count(status)>0)  // if this status is auto-advanced to another status
                {
                    std::string nextStatus = statusSequences.at(status).next;
                    std::tuple<double, double> statusDuration= statusSequences.at(nextStatus).duration;
                    set_status(*it, t, nextStatus, statusDuration);
                }


                std::iter_swap(lo_it, it); // switch expired farm into low position
                lo_it++; // shift lo placemarker forward
                inStatusList.lo++; // update stored integer to match lo_it position (as integer)
            }

        }
    } // end "if there are farms in list that haven't expired"

}

/// Overloaded for Region_status. Unlike with premises, regions are reported via the
/// report_countyAndState function.
void Status_manager::update(int t, std::string status, statusList<Region_status*>& inStatusList)
{
	// check from iterator "lo" to end, switch expired farms to before "lo"
	if (inStatusList.lo < inStatusList.units.size()){ // not all farms have expired
		std::vector<Region_status*>::iterator lo_it = inStatusList.units.begin();
		std::advance(lo_it, inStatusList.lo); // move iterator to lo

		for (auto it = lo_it; it < inStatusList.units.end(); it++){ // check each farm* between lo and end
				if ( (*it)->Region_status::get_end(status) == t ){ // if validity of this status expires for this farm today
				// set to next status
				if (statusSequences.count(status)>0){ // if this status is auto-advanced to another status
					std::string nextStatus = statusSequences.at(status).next;
					std::tuple<double, double> duration = statusSequences.at(nextStatus).duration;
					std::string id = (*it)->Region::get_id();
					std::string regionType = (*it)->Region::get_type();
					set_regionStatus(id, regionType, t, nextStatus, duration);
				}
				std::iter_swap(lo_it, it); // switch expired region into low position
				lo_it++; // shift lo placemarker forward
				inStatusList.lo++; // update stored integer to match lo_it position (as integer)
				}
		}
	} // end "if there are farms in list that haven't expired"

}



/// Update disease statuses (wrapper for update() function)
/// \param[in] t Timestep
void Status_manager::updateDisease(int t)
{
	// 'sus'->'exp' transition determined by local spread, shipping, +control
	if (diseaseStatuses.count("exp")>0){
		update(t, "exp", diseaseStatuses.at("exp")); // any 'exp' expiring -> 'inf'
	}
	if (diseaseStatuses.count("inf")>0){
	update(t, "inf", diseaseStatuses.at("inf")); // any 'inf' expiring -> 'imm'
	}
	// 'imm' is a permanent status
}


/// Update file and diagnostics statuses (wrapper for update() function)
/// \param[in] t Timestep
void Status_manager::updateFile(int t)
{
	// 'sus'->'exposed' transition determined by local spread, shipping, +control
	if (fileStatuses.count("exposed")>0){
		update(t, "exposed", fileStatuses.at("exposed")); // any 'exposed' expiring -> reported+waitlisted
	}
	// 'reported' is a permanent status, no updates performe
}


/// Check time for time dependent resource boost (i.e. for additional vaccine doses)
void Status_manager::update_ControlResources(int t)
{
		for (auto& ct:(parameters->controlTypes)){
		std::string c_type = ct;
 	   	std::string c_constraint = (parameters->constraintFunctions).at(ct);
 	   	
		if((c_type.compare("vax")==0) && (c_constraint.compare("nationalLimit")==0)){
			Control_resource* cr = controlResources->at("vax").at("all");
				if (controlResourceLevels.count(cr)==0){ //if levels for this control resource are not yet being tracked
					controlResourceLevels[cr] = cr->Control_resource::get_capacity();
				}else{

			if(firstRepTime > 0){

				int timeToCheck = t-firstRepTime; //Current time minus time of first report

				if(controlReleaseSchedule->at("vax").count(timeToCheck)>0){
					int boost = controlReleaseSchedule->at("vax").at(timeToCheck);
					controlResourceLevels[cr] += boost;
					if(verbose>1){
						std::cout << "controlResourceLevels "<< controlResourceLevels[cr]<<std::endl;
					}
					if(verbose>1){
						std::cout << "SM: Boost of " <<boost<<" at "<<timeToCheck<< std::endl;
					}
				} }else{if(verbose>1)
					std::cout << "SM: First reported time not yet recorded" << std::endl;
					}
			//	std::cout<<"Available resources at this time point are for "<<c_type<<" are "<<cr->Control_resource::get_capacity()<<std::endl;
		} 
	}
	}
}

/// Update file and control statuses (wrapper for update() function)
/// \param[in] t Timestep
void Status_manager::updateControl(int t)
{
	// update control statuses
	for (auto& ct:*allControlTypes){ // should be all control types, including regional
		// 'waitlisted' & 'implemented' transitions determined by controlRules, constraints
		if (controlStatuses.count("implemented."+ct.first)>0){
			update(t, "implemented."+ct.first, controlStatuses.at("implemented."+ct.first)); // any 'implemented' expiring -> 'effective'
		} else if (regionControlStatuses.count("implemented."+ct.first)>0){
			update(t, "implemented."+ct.first, regionControlStatuses.at("implemented."+ct.first)); // any 'implemented' expiring -> 'effective'
		}
		if (controlStatuses.count("effective."+ct.first)>0){
			update(t, "effective."+ct.first, controlStatuses.at("effective."+ct.first)); // any 'effective' expiring -> 'inactive'
		} else if (regionControlStatuses.count("effective."+ct.first)>0){
			update(t, "effective."+ct.first, regionControlStatuses.at("effective."+ct.first)); // any 'effective' expiring -> 'inactive'
		}
		// 'inactive' is a permanent status
	}

}


/// Update file and diagnostics statuses (wrapper for update() function)
/// \param[in] t Timestep
void Status_manager::updateDiagnostics(int t)
{
	// update diagnostic statuses
	for (auto& dt:*allDiagnosticTypes){ // should be all diagnostic types
		if (diagnosticStatuses.count("started."+dt.first)>0){
			update(t, "started."+dt.first, diagnosticStatuses.at("started."+dt.first)); // any 'started' expiring -> 'complete'
		}
		if (diagnosticStatuses.count("complete."+dt.first)>0){
			update(t, "complete."+dt.first, diagnosticStatuses.at("complete."+dt.first)); // any 'started' expiring -> 'complete'
		}
	}

}

///Evaluate the diagnostic test and save the positive or negative result to Prem_Status
/// \parm[in] t Timestep
void Status_manager::eval_testResult(int t, int quarter_idx){

	for (auto& dt:toTest){ // should be all diagnostic types
		if(dt.second.size()>0){
		std::string d_type = dt.first;
		std::vector<Farm*> toTest;


		for (auto& f:dt.second){
			bool testPositive = eval_probHerdPositive(f, t, quarter_idx, d_type);
			int fid = f->Farm::get_id();

			changedStatus.at(fid)->Prem_status::add_testResult(d_type, 0);

			if(testPositive == 0){
				if (verbose>1){std::cout<<"Diagnostic test: negative"<<std::endl;}
				changedStatus.at(fid)->Prem_status::add_testResult(d_type, 0);
			}else if(testPositive == 1){
				if (verbose>1){std::cout<<"Diagnostic test: positive"<<std::endl;}
				changedStatus.at(fid)->Prem_status::add_testResult(d_type, 1);

			}
		}
	}
  }
}

/// Stochastically determines if the diagnostic test is positive or negative at a given premises.
/// Returns 'true' if test is positive, 'false' if it is not.
/// \param[in] farm id
bool Status_manager::eval_probHerdPositive(Farm* f, int t, int quarter_idx, std::string d_type){
	bool output = false; // by default assume test is negative
    int fid = f->Farm::get_id();
    if (changedStatus.count(fid)>0)//Diagnostics begun on focal farm but not complete
    {

		std::tuple<double, double> sens = allDiagnosticTypes->at(d_type)->sensitivity;

			if(parameters->infectionType == InfectionType::FMD &&
               (changedStatus.at(fid)->Prem_status::get_diseaseStatus()=="inf" ||
                changedStatus.at(fid)->Prem_status::get_diseaseStatus()=="exp" ||
				changedStatus.at(fid)->Prem_status::get_diseaseStatus()=="imm"))
            {
                double n = 0;
                if(parameters->partial==1){
                     n = changedStatus.at(fid)->Prem_status::get_inf_partial_as_vaccinated(t, quarter_idx, parameters, this->get_normInf_map()); //Returns the number of unvaccinated infectious animals. If vaccine hasn't been deployed, the number of unvaccinated animals equals total herd size.
                } else if (parameters->partial==0){
                     n = changedStatus.at(fid)->Farm::get_size_allSpecies();
                }

				double normDraw = draw_norm(std::get<0>(sens), std::get<1>(sens));
				double bracket =1-normDraw;
				double probNCatNeg = pow(bracket, n);
				double probHerdPositive = 1-probNCatNeg;
		 		double random = uniform_rand();
		 		if (random < probHerdPositive)
                {
		 			output = true;
		 		}
		 	}
		 	else if(parameters->infectionType == InfectionType::BTB &&
                    changedStatus.at(fid)->Prem_status::get_diseaseStatus()=="btb")
            {
                int n_infected = changedStatus.at(fid)->btb_get_N_infected();
                double normDraw = draw_norm(std::get<0>(sens), std::get<1>(sens));
				int n_positive = draw_binom(n_infected, normDraw);
		 		if (n_positive > 0)
                {
		 			output = true;
		 		}
            }
            else //diseaseStatus == "sus"
            {
                if (verbose>1){std::cout<<"SM: eval_probHerdPositive no test evaluation performed, farm is susceptible"<<std::endl;}
                output = false;
            }

	}
	return output;
}

void Status_manager::eval_slaughterSurveillance(int t, int quarter_idx, int replicate, std::vector<Farm*>prems)
{
    if(!btb_slaughter_ofs.is_open())
    {
        btb_slaughter_ofs.open(btb_slaughter_ofn, std::ios::app);
        if(!btb_slaughter_ofs){std::cout << "Failed to open " << btb_slaughter_ofn << std::endl; exit(EXIT_FAILURE);}
        //Write the header.
        if(replicate == 1)
        {
            btb_slaughter_ofs << "replicate" << "\t" << "timestep" << "\t" << "oPremId" << "\t" << "oCommodity" << "\t" << "oPremClass" << "\t"
                              << "oPremSize" << "\t" << "oPremNInf" << "\t" << "nToSlaughter" << "\t" << "nInfToSlaughter" << "\t"
                              << "nDetectableLesions" << "\t" << "nInfAnimalsDetected" << "\t" << "tracebackSuccessful" << "\t" << "facilityId" << "\t"
                              << "infectionTracedTo" << std::endl;
        }
    }
    for(Farm* f : prems)
    {
        Prem_status* pst = changedStatus.at(f->get_id());
        if(pst->get_fileStatus() == "reported")
        {
            //Reported prems don't sent to slaughter.
            continue;
        }
        double slaughter_shipment_factor = f->get_slaughter_shipment_factor();
        if(slaughter_shipment_factor < 0.0)
        {
            //Slaughter shipment factor hasn't been set yet and has the default value of -1.0.
            slaughter_shipment_factor = gridManager->get_slaughter_shipment_factor(f->get_farm_type()->get_species(),
                                                                                   f->get_prem_class()->tag,
                                                                                   f->get_USAMMv3_unbinned_size());
            f->set_slaughter_shipment_factor(slaughter_shipment_factor);
        }

        double slaughter_shipment_rate = f->get_shipment_rate() * slaughter_shipment_factor; //This is the monthly rate for this quarter.
        if(slaughter_shipment_rate > 0.0)
        {
            int n_slaughter_shipments = draw_poisson(slaughter_shipment_rate); //Generate the number of slaughter shipmetns that occur this time step.
            if(n_slaughter_shipments > 0)
            {
                double slaughter_fraction = f->get_yearly_slaughter_fraction();
                if(slaughter_fraction < 0.0)
                {
                    //Slaughter fraction hasn't been set yet and has the default value of -1.0.
                    slaughter_fraction = gridManager->get_slaughter_fraction(f->get_farm_type()->get_species(),
                                                                             f->get_prem_class()->tag);
                    f->set_yearly_slaughter_fraction(slaughter_fraction);
                }

                double yearly_sl_rate = slaughter_shipment_rate * 12; //This is also the expected number of yearly slaughter shipments given this quarter's monthly rate.
                double fraction_per_sl_shipment = slaughter_fraction / yearly_sl_rate; //This is the expected proportion of the inventory that is sent to slaughter on each slaughter shipment given the quarter's shipping rate.


                Prem_status* farm_triggered = nullptr;
                for(int sl_shp_idx = 0; sl_shp_idx<n_slaughter_shipments; ++sl_shp_idx)
                {
                    int n_to_slaughter = 0;
                    int n_infected = 0;
                    int n_detectable_lesions = 0;
                    int confirmed_inf_animals = 0;
                    int traceback_successful = 0;
                    int facility_id = -1;
                    if(fraction_per_sl_shipment >= 1.0)
                    {
                        n_to_slaughter = f->get_size_allSpecies();
                    }
                    else
                    {
                        int premsize = f->get_size_allSpecies();
                        if(premsize == 1)
                        {
                            n_to_slaughter = 1;
                        }
                        else
                        {
                            n_to_slaughter = 1 + draw_binom(f->get_size_allSpecies()-1, fraction_per_sl_shipment); //To prevent sending zero animals.
                        }

                    }

                    n_infected = pst->btb_subtract_slaughter_shipment(n_to_slaughter); //Hypergeometric sample from the different infection stages.
                    if(n_infected > 0)
                    {
                        //Determine how many lesions are detectable. Same probability regardless of facility.
                        n_detectable_lesions = 0;
                        for(int animal_i=0; animal_i<n_infected; ++animal_i)
                        {
                            double p_detectable = draw_beta(30, 10);
                            if(uniform_rand() < p_detectable)
                            {
                                ++n_detectable_lesions;
                            }
                        }

                        //Determine how many of the lesions are detected and submitted.
                        if(n_detectable_lesions > 0)
                        {
                            //Get the destination slaughter facility.
                            County* o_county = f->get_parent_county();
                            if(!o_county->is_set_slaughter_probs())
                            {
                                int o_county_fips = o_county->get_fips_code();
                                o_county->set_slaughter_probs(gridManager->get_slaughtershed_probabilities(o_county_fips),
                                                              gridManager->get_slaughter_facility_ids());
                            }
                            facility_id = o_county->get_slaughter_destination();
                            if(facility_id != -1)
                            {
                                //Now determine slaughter facility/chance of submitting a lesion and then if that lesion
                                //tests positive and if so a traceback is successful.
                                const std::vector<double>& alpha_beta = gridManager->get_slaughter_facility_alpha_beta(facility_id);
                                confirmed_inf_animals = 0;
                                for(int lesion_i=0; lesion_i<n_detectable_lesions; ++lesion_i)
                                {
                                    double p_detect = draw_beta(alpha_beta[0], alpha_beta[1]);//Draw one random probability for each animal checked.
                                    if(uniform_rand() < p_detect)
                                    {
                                        //Lesion detected and submitted.
                                        double p_test_positive = draw_beta(46, 5);
                                        if(uniform_rand() < p_test_positive)
                                        {
                                            ++confirmed_inf_animals;
                                        }
                                    }
                                }

                                if(confirmed_inf_animals > 0)
                                {
                                    //At least one lesion tested positive for bTB, try to trace it back to the prem.
                                    double p_traceback = draw_beta(45, 51);
                                    if(uniform_rand() < p_traceback)
                                    {
                                        traceback_successful = 1;
                                        if(pst->get_prem_class()->tag == "Frm")
                                        {
                                            //Infected prem found, trigger diagnostics measures.
                                            farm_triggered = pst;
                                        }
                                        else
                                        {
                                            //It's not a farm, so look through the received shipments to see if
                                            //a farm can be found.
                                            farm_triggered = btb_make_trace(pst); //The prem was reported so make a trace in order to find the source.
                                        }
                                    }
                                }
                            }
                        }
                    }

                    int farm_triggered_id = -1;
                    if(farm_triggered != nullptr and farm_triggered->get_prem_class()->tag == "Frm")
                    {
                        farm_triggered_id = farm_triggered->get_id();
                    }

                    //Print the slaughter shipment to file.
                    std::string commodity = f->get_farm_type()->get_species();
                    btb_slaughter_ofs << replicate << "\t"
                                      << t << "\t"
                                      << f->get_id() << "\t"
                                      << commodity << "\t"
                                      << f->get_prem_class()->tag << "\t"
                                      << f->get_USAMMv3_unbinned_size() << "\t"
                                      << pst->btb_get_N_infected() << "\t"
                                      << n_to_slaughter << "\t"
                                      << n_infected << "\t"
                                      << n_detectable_lesions << "\t"
                                      << confirmed_inf_animals << "\t"
                                      << traceback_successful << "\t"
                                      << facility_id << "\t"
                                      << farm_triggered_id << std::endl;
                    if(farm_triggered != nullptr and farm_triggered->get_prem_class()->tag == "Frm")
                    {
                        if(parameters->testSuspects_on)
                        {
                            set_status(farm_triggered, t, "suspected", std::make_tuple<double, double>(0.0, 0.0));
                        }
                        else
                        {
                            set_status(farm_triggered, t, "reported", pastEndTime);
                        }

                        break; //No need to test the rest of the slaughter shipments, diagnostics already triggered for this prem.
                    }
                }
            }
        }
    }
}

void Status_manager::wipe_markets()
{
    for(auto& id_pst_pair : changedStatus)
    {
        Prem_status* pst = id_pst_pair.second;
        if(pst->get_prem_class()->tag == "Mkt" and pst->get_diseaseStatus() == "btb")
        {
            pst->btb_remove_all_infected();
        }
    }
}

///Stochastically determines if the vaccinated animals that are present on one or both
///of the premises suppress the transmission event. If the suppression is successful
///and transmission is prevented whichever of the boolean reference paramters that is
///appropriate will be set to false: if the origin herd has been vaccinated, transmission
///is prevented from that herd and transPrevented is true; if the receiveing herd is
///vaccinated, exposure is prevented and expPrevented is true; if both are vaccinated
///there is no easy way to determine which vaccinateion event prevented exposure and
///both transPrevented and expPrevented are set to true.
void Status_manager::eval_premExpPrevByVaccination(Farm* ofarm, Farm* dfarm, double trueP, int t,
                                                   int quarter_idx, bool& transPrevented,
                                                   bool& expPrevented)
{
    int ofid = ofarm->get_id();
    int dfid = dfarm->get_id();

    //First determine which of the farms that are vaccinated.
    bool ofarm_isvaxed = false;
    Prem_status* opst = nullptr;
    try{
        opst = changedStatus.at(ofid);
        ofarm_isvaxed = opst->get_isVaccinated();
    }
    catch(...) {} //Just continue if there is no vaccine in effect.

    //Same as above for destination farm.
    bool dfarm_isvaxed = false;
    Prem_status* dpst = nullptr;
    try{
        dpst = changedStatus.at(dfid);
        dfarm_isvaxed = dpst->get_isVaccinated();
    }
    catch(...) {}

    if(!ofarm_isvaxed and !dfarm_isvaxed)
    {
        //No farm is vaccinated so skip rest
        transPrevented = false;
        expPrevented = false;
        return;
    }
    //Now we calculate the probability of transmission with completely unvaccinated herds
    //and the probability with vaccinated herds. Transmission has already been detemined
    //to take place under the condition that there is no control in place so now we need
    //the conditional probability of infection happening to the vaccinated herds given
    //that transmission would have taken place if they were both unvaccinated:
    // p_actual = p_vax / p_novax. For the prob that the transmission is blocked
    // we take 1-p_actual = p_blocked. For efficiency the true probability of transmission
    //that caused the infection, based on unvaccinated populations iss deconstructed into
    //the riginal rate and that rate is adjusted to use the vaccinated populations instead
    //and then converted back into a (lower) probability.
    //The rate for unvaccinated herds is r_ij = S*T*Ni*Nj*K and the rate for
    //vaccinated populations is rv_ij = S*T*Nvi*Nvj*K, so multiply r_ij by the ratios
    //Nvi / Ni and Nvj / Nj to get rv_ij.

    //Determine what to be used for infectiousness of the origin farm.
    double ofarm_inf_unvaxed; //Transmissibility of origin farm.
    double ofarm_inf_vaxed; //Transmissibility of origin farm using the susceptible parts of the herd when vaccination is in place.
    double dfarm_sus_unvaxed = dfarm->get_sus(quarter_idx); //Same but susceptibility of destination farm. Will always be based on full pop regardless of partial transmission or not.
    double dfarm_sus_vaxed = dfarm_sus_unvaxed;
    auto& normInf_map = this->get_normInf_map();
    auto& normSus_map = this->get_normSus_map();
    if(parameters->partial)
    {
        ofarm_inf_unvaxed = opst->get_inf_partial_as_unvaccinated(t, quarter_idx, parameters, normInf_map);
        if(ofarm_isvaxed)
        {
            ofarm_inf_vaxed = opst->get_inf_partial_as_vaccinated(t, quarter_idx, parameters, normInf_map);
        }
        else
        {
            ofarm_inf_vaxed = ofarm_inf_unvaxed;
        }
    }
    else
    {
        ofarm_inf_unvaxed = opst->get_inf(quarter_idx);
        ofarm_inf_vaxed = 0.0;
        for(const auto& sp_N_pair : opst->get_currentSizeUnvaccinated())
        {
            int N = sp_N_pair.second;
            if(N > 0)
            {
                ofarm_inf_vaxed += normInf_map.at(sp_N_pair.first).at(quarter_idx)  * std::pow(N, parameters->infExponents.at(sp_N_pair.first));
            }
        }
    }

    if(dfarm_isvaxed)
    {
        dfarm_sus_vaxed = 0.0;
        for(const auto& sp_N_pair : dpst->get_currentSizeUnvaccinated())
        {
            int N = sp_N_pair.second;
            if(N > 0)
            {
                dfarm_sus_vaxed += normSus_map.at(sp_N_pair.first).at(quarter_idx)  * std::pow(N, parameters->susExponents.at(sp_N_pair.first));
            }
        }
    }

    double original_rate = -std::log(-(trueP-1));
    double vaccinated_rate = original_rate * ( (ofarm_inf_vaxed * dfarm_sus_vaxed) / (ofarm_inf_unvaxed * dfarm_sus_unvaxed) );
    double p_exp_vaccinated = oneMinusExp(-vaccinated_rate * parameters->days_per_timestep); //This is the probability that the transmission event would have occurred with vaccinated animals.
    double p_conditional = p_exp_vaccinated / trueP; //This is the probability that the transmission takes place between vaccinated herds given that it would have taken place with unvaccinated herds.
    double p_prevent = 1 - p_conditional; //This is the complement to the above (that the transmission does not take place between vaccinated herds given...).

    //Now to se if the transmission is blocked
    double rval = uniform_rand();
    if (rval <= p_prevent) //Exposure blocked.
    {
        if(ofarm_isvaxed and !dfarm_isvaxed)
        {
            transPrevented = true;
        }
        else if(!ofarm_isvaxed and dfarm_isvaxed)
        {
            expPrevented = true;
        }
        else if(ofarm_isvaxed and dfarm_isvaxed)
        {
            transPrevented = true;
            expPrevented = true;
        }
    }
    else
    {
        //Else the transmission wasn't blocked by the precense of vaccinated animals
        //on one or both of the farms and both bools are true, both transmission from
        //the infectious farm and exposure at the receiveing farm takes place.
        transPrevented = false;
        expPrevented = false;
    }
}


/// Stochastically determines if the control effective at a given premises actually
/// prevents a transmission event. Returns 'true' if transmission occurs, 'false' if it does not.
/// Used with evaluating local spread and shipments escaping ban.
/// param[in] f Farm* to check
bool Status_manager::eval_premTransmission(Farm* f){
	bool output = 1; // by default assume transmission occurs
	int fid = f->Farm::get_id();
	if (changedStatus.count(fid)>0){
		double pNoTransmission = changedStatus.at(fid)->Prem_status::get_probPreventTransmission();
		double random = uniform_rand();
		if (random < pNoTransmission){  // farm does not become exposed due to control
			output = 0;
		}
	} // otherwise, farm has not had any status change, including control, so will be exposed
	return output;
}

/// Stochastically determines if the control effective at a given premises actually
/// prevents an exposure event. Returns 'true' if exposure occurs, 'false' if it does not.
/// Used with evaluating local spread and shipments escaping ban.
/// param[in] f Farm* to check
bool Status_manager::eval_premExposure(Farm* f){

	bool output = 1; // by default assume exposure occurs
	int fid = f->Farm::get_id();
	if (changedStatus.count(fid)>0){
		double pNotExposed = changedStatus.at(fid)->Prem_status::get_probPreventExposure();
		double random = uniform_rand();
		if (random < pNotExposed){  // farm does not become exposed due to control probability
			output = 0;
		}
	} // otherwise, farm has not had any status change, including control, so will be exposed
	return output;
}

/// Evaluates all premises-level control measures impacting local spread from locally-exposed and
/// shipping-exposed farms, and exposes those farms where exposure is not prevented. If
/// any farms are exposed through multiple local or shipment events, it is
/// evaluated for each exposure independently. If
/// exposure is not prevented in at least one of those evaluations, the farm becomes exposed.
/// \param[in] t Current timestep at which exposure/blocking events are reported
void Status_manager::eval_exposure(int t, int quarter_idx)
// Exposures with sources of infection stored in expForEval
{

	int localCount = 0;
	int shipCount = 0;
	std::vector<std::pair<Farm*, int>> toExpose;
	for (auto& e:exposureForEval){ // e = tuple of destination Farm, source Farm, route, true probability of infection
		Farm* destination = std::get<0>(e);
		Farm* origin = std::get<1>(e);
		int route = std::get<2>(e);
		double expP = std::get<3>(e);
		std::string prevented = "none"; // default value - no prevention, exposure occurs
        int latency_for_this_exposure = -1; //Default -1 means use latency from the config file.
		if (parameters->control_on == 0){
			add_premSource(t, destination, origin, route, prevented);
			if(route == 1) //Exposure through shipment
            {
                Prem_status* origin_prem_status = changedStatus.at(origin->Farm::get_id());
                if(origin_prem_status->get_diseaseStatus().compare("exp") == 0) //The origin farm is exposed and not infectious yet. The destination farm will inherit the origins' latency to infectiousness.
                {
                    latency_for_this_exposure = std::max(1, origin_prem_status->get_end("exp") - t); //Latency is inherited. Must never be 0 as the "exp" status is skipped over with bad consequences in the set_status if that's the case.
                }
                else if(origin_prem_status->get_diseaseStatus().compare("inf") == 0) //Not sure if this should be a feature. Need to decide if bringing infectious animals onto a farm leads to instant infectiousness. /Stefan
                {
                    latency_for_this_exposure = 1; //Instant infectiousness.
                }
            }
            std::pair<Farm*, int> exp_pair = std::make_pair(destination, latency_for_this_exposure);
			toExpose.emplace_back(exp_pair);
			if (route==0){++localCount;}
			else if (route==1){++shipCount;}

		} else if (parameters->control_on == 1){
			// determine if origin farm is controlled and transmits

			//First check if vaccination present in one or both herds prevents transmission/exposure.
			//If one and not the other herd is vaccinated we can determine which of the herds that prevented
			//the transmission or exposure (depending on which of the two herds were vaccinated). If both herds
			//are vaccinated it's not possible to determine the probability of one or the other actually contributing
			//to the preventing of the infection. This is because the imperfect vaccinaton is applied to
			//individual animals and therefore affects the rate of infection which is a function of both herd sizes
			//and not just one of them. Therefore if both are vaccinated they both contribute to lowering the rate
			//compared to what it would have been if they were unvaccinated. Look at it this way: if one herd consisted
			//only two animals, and letting the susceptibility and transmissibility
			//constants S=T=0.1 (and assuming a linear relationship between herd size and susc. & transm.) the
			//probability of infecting another herd of say 10 unvaccinated animals would be p = 1-e^-(2*10*0.01) = 0.181.
			//Lets say one of the two animals are vaccinated, then p = = 1-e^-(1*10*0.01) = 0.095, a reduction by 47.5%.
			//But what if the second herd consisted of 1000 animals, then the same reduction would only be 26.9 %. Clearly
			//the probability of vaccine blocking transmission is dependent of both premises' vaccnated herd size.
			//For this reason, if both premises have vaccinated animals and the reduced rate leads to infection being
			//blocked we let the outcome be that both transmission and exposure was blocked.

			bool localTransmissionOccurs = eval_premTransmission(origin);
			bool localExposureOccurs = eval_premExposure(destination);

			bool vaxPreventsTransmission = false;
            bool vaxPreventsExposure = false;
			if(localTransmissionOccurs or localExposureOccurs) //No use if event is prevented by other means.
            {
                eval_premExpPrevByVaccination(origin, destination, expP, t, quarter_idx,
                                              vaxPreventsTransmission, vaxPreventsExposure);
            }

            localTransmissionOccurs = localTransmissionOccurs && !vaxPreventsTransmission; //Multiplying bools is equivalent to multiplying ones and zeros (same as AND).
            localExposureOccurs = localExposureOccurs && !vaxPreventsExposure;

			if (!localTransmissionOccurs){
				if (verbose>1){std::cout<<"Transmission blocked."<<std::endl;}
				prevented = "src:"+vecToCommaSepString(changedStatus.at(origin->Farm::get_id())->Prem_status::get_controlStatuses()); // determining which of multiple control types could get complicated
				add_premSource(t, destination, origin, route, prevented);
			}

			if (!localExposureOccurs){
				if (verbose>1){std::cout<<"Exposure blocked."<<std::endl;}
				prevented = "exp:"+vecToCommaSepString(changedStatus.at(destination->Farm::get_id())->Prem_status::get_controlStatuses()); // determining which of multiple control types could get complicated
				add_premSource(t, destination, origin, route, prevented);
			}
			else if (localTransmissionOccurs && localExposureOccurs){
				add_premSource(t, destination, origin, route, prevented);
				if(route == 1) //Exposure through shipment
                {
                    Prem_status* origin_prem_status = changedStatus.at(origin->Farm::get_id());
                    if(origin_prem_status->get_diseaseStatus().compare("exp") == 0) //The origin farm is exposed and not infectious yet. The destination farm will inherit the origins' latency to infectiousness.
                    {
                        latency_for_this_exposure = origin_prem_status->get_end("exp") - t; //Latency is inherited.
                    }
                    else if(origin_prem_status->get_diseaseStatus().compare("inf") == 0) //Not sure if this should be a feature. Need to decide if bringing infectious animals onto a farm leads to instant infectiousness. /Stefan
                    {
                      latency_for_this_exposure = 0; //Instant infectiousness.
                    }
                }
                std::pair<Farm*, int> exp_pair = std::make_pair(destination, latency_for_this_exposure);
				toExpose.emplace_back(exp_pair);
				if (route==0){++localCount;}
				if (route==1){++shipCount;}
			}
		} // end "if control on"
	} // end "for each exposure to evaluate"
	if(verbose>0){
		std::cout<<localCount<<" local infections (may include double-counted)"<<std::endl;
		std::cout<<shipCount<<" shipping infections (may include local)"<<std::endl;
	}
	expose(toExpose, t);
	exposureForEval.clear();
}

/// Filters exposures that are prevented via control (shipBan), records sources of
/// exposure in Prem_status, returns farms that would be exposed. Output may include farms
/// that are already infected. Premises-level control is checked at a later step
/// (Status_manager::eval_exposure())
void Status_manager::filter_shipments(std::vector<Shipment*>& ships, int time)
// filters for any control measures impacting shipment spread
// records sources of infection as shipment
// shipment has t, farm origID, farm destID, origin FIPS, dest FIPS, species, ban
{
	for (auto& s:ships){
	// fill in fields t, origin, destination
		s->timestep = time; // set time of shipment
		Farm* destination = s->dPrem;

		// Only evaluate exposure and control if destination is susceptible
		if (getAny_diseaseStatus(destination).compare("sus")==0 or
            getAny_diseaseStatus(destination).compare("btb")==0){
			Farm* origin = s->oPrem;
			bool exposeDestination = true; // default assumption, control will turn this off
			if(parameters->market_shipment_exposure_prob < 1.0 and origin->get_prem_class()->tag == "Mkt")
            {
                if(uniform_rand() > parameters->market_shipment_exposure_prob)
                {
                    //This shipment wasn't infected. Just skip to next.
                    continue;
                }
            }
			// Exposure does NOT happen if shipping bans are effective and realized:
			if (parameters->control_on == true && allControlTypes->count("shipBan")>0){
				std::string county_fips = origin->Farm::get_parent_county()->Region::get_id();
				std::string state = origin->Farm::get_parent_state()->Region::get_id();
				std::string shipBanScale = allControlTypes->at("shipBan")->scale;
				if (shipBanScale.compare("county")==0 && changedCoStatus.count(county_fips) > 0){
					double pBan = changedCoStatus.at(county_fips)->Region_status::get_probPreventExposure();
					if (pBan > 0){ // county shipBan is effective
						double random = uniform_rand();
						if (random < pBan){ // county shipBan is realized
							s->ban = 1;
							exposeDestination = false;
if(verbose>1){std::cout<<"SM::filter_shipments: Shipment prevented by ban"<<std::endl;}
						}
					} // end "if county shipBan is effective"
				} else if (shipBanScale.compare("state")==0 && changedStateStatus.count(state)>0){
					double pBan = changedStateStatus.at(state)->Region_status::get_probPreventExposure();
					if (pBan>0){ // state shipBan is effective
						double random = uniform_rand();
						if (random < pBan){ // state shipBan is realized
							s->ban = 1;
							exposeDestination = false;
if(verbose>1){std::cout<<"SM::filter_shipments: Shipment prevented by ban"<<std::endl;}
						} else {
						if(verbose>1){std::cout<<"Shipment avoided ban"<<std::endl;}}
					} // end "if state shipBan is effective"
				}	// end "if county or state level ban"
				// No other regional shipment bans currently implemented
			} // end "if control_on"

			// record exposure source in Prem_status
			if (changedStatus.count(origin->Farm::get_id())<1){
				std::cout<<"ERROR: In Status_manager::filter_shipments: Assumed infectious shipment originated from farm "<<origin->Farm::get_id()<<", but farm not recorded as infectious. Exiting...";
				exit(EXIT_FAILURE);
			}

			if (!exposeDestination){ // exposure was blocked by shipBan, record, but nothing further will happen with this premises
				add_premSource(time, destination, changedStatus.at(origin->Farm::get_id()), 1, "shipBan");
			} else if (exposeDestination){ // control not on or not realized
				// evaluate for prem-level control, exposure
				add_premForEval(destination, origin, 1, 0.0, s); // stores info temporarily in Prem_status
			}

		} // end "if destination is susceptible"
	} // end "for each shipment"
}

/// Sets farm diseaseStatus to "exposed" if status is currently "susceptible" (any
/// subsequent attempts to expose a farm will be ignored). If control is on, also
/// sets up farm fileStatus for reporting, based on dangerousContact status.
/// the input vector consists of pairs of farm to be exposed and explicit latency
/// for that exposure, if that latency is negative, the default will be used (as defined
/// in the config file), if it is set to zero the farm will be instantly infectious. This
/// is to enable exposure by exposed shipments and instant infectiousness by moving
/// infectious animals onto a farm.
void Status_manager::expose(std::vector<std::pair<Farm*, int>>& farms, int t)
{
	for (auto& farm_latency_pair : farms){
        Farm* f = farm_latency_pair.first; //Farm to which to apply status change.
        int explicit_latency = farm_latency_pair.second; //Explicit latency for this exposure.
		// check if farm is susceptible
		if (getAny_diseaseStatus(f).compare("sus")==0){
		    std::tuple<double, double> latency_parameters = parameters->latencyParams; //Default parameters are those specified in config file
		    if(explicit_latency > -1) //A specific latency will be used instead of the default.
            {
                std::get<0>(latency_parameters) = explicit_latency; //Set mean to explicit latency...
                std::get<1>(latency_parameters) = 0.0; //...and variance to zero.
            }
			set_status(f,t,"exp",latency_parameters); // disease status
			if (parameters->control_on == true || parameters->diagnostic_on == true){
			// determine time to reporting based on dangerousContact status
			  if (parameters->testSuspects_on == false){
				if (getAny_fileStatus(f).compare("dangerousContact") == 0){
					set_status(f,t,"exposed",parameters->dcReportLag); // altered report time for DCs
				} else if (getAny_fileStatus(f).compare("notDangerousContact") == 0){
					set_status(f,t,"exposed",parameters->nonDCReportLag);
				}
			 } else if (parameters->testSuspects_on == true){
			// determine time to reporting based on dangerousContact status
				if (getAny_fileStatus(f).compare("notDangerousContact") == 0){
					set_status(f,t,"exposed",parameters->nonDCSuspectLag);
				} else if (getAny_fileStatus(f).compare("dangerousContact") == 0){
					set_status(f,t,"exposed",parameters->dcSuspectLag);
				}
			}//end "if test suspects on"
		}// end "if control on"
		} // end "if farm is susceptible"
	} // end "for each farm"
}



/// For a vector of given farms, reports the parent county and state.
void Status_manager::report_countyAndState(Farm* f, int t)
{
	County* parentCounty = f->Farm::get_parent_county();
	std::string parentFips = parentCounty->Region::get_id();
	State* parentState = f->Farm::get_parent_state();
	std::string stateName = parentState->Region::get_id();

	if (changedCoStatus.count(parentFips) == 0){
		changedCoStatus[parentFips] = new Region_status(parentCounty);
	}
	if (changedCoStatus.at(parentFips)->Region_status::is_reported() == false){
		changedCoStatus.at(parentFips)->Region_status::report();
		reportedCounties.emplace_back(parentFips);
		newCoReports.emplace_back(changedCoStatus.at(parentFips));
	}

	if (changedStateStatus.count(stateName) == 0){
		changedStateStatus[stateName] = new Region_status(parentState);
	}
	if (changedStateStatus.at(stateName)->Region_status::is_reported() == false){
		changedStateStatus.at(stateName)->Region_status::report();
		reportedStates.emplace_back(stateName);
		newStateReports.emplace_back(changedStateStatus.at(stateName));
	}
}


/// Passes conditions about current status of outbreak (number of premises infected, etc)
/// and checks if control rules indicate a control action should be triggered. Adds
/// premises or regions to appropriate waitlists.
void Status_manager::add_waitlistMembers(int t)
{
if (verbose >1 && newPremReports.size()+newCoReports.size()+newStateReports.size()>0){
std::cout<<"SM::add_waitlistMembers: Checking for control triggers: "<<newPremReports.size()<<" new prem reports, "<<
newCoReports.size()<<" new county reports, "<< newStateReports.size()<<" new state reports"<<std::endl;
}
	// statsForControlTriggers and waitlistGroup structs defined in Control_manager
	statsForControlRules stats{&newPremReports, &newCoReports, &newStateReports, percentIncrease, decrease, controlsTriggered, controlResourceLevels, waitlistFlag, t};
	std::vector<waitlistGroup> instructions;
	controlManager->check_controlRules(stats, instructions); // writes waitlist info to instructions
	// add to waitlist according to instructions
	for (auto& i:instructions){

		std::string controlType = i.controlType;
		std::string controlScale = i.scaleType;
		double controlSize = i.premList.size();

		//start tracking what controls have been triggered
		if(controlSize>0 && (std::find(controlsTriggered.begin(), controlsTriggered.end(), controlType) == controlsTriggered.end())){		
			controlsTriggered.emplace_back(controlType);
		}

//std::cout<<"replaceControl check. Is it on? "<<i.replaceControl<<std::endl;
//precursor to replace control actions (if cull has been called already this would allow us replace cull with a different type of cull and turn off the old type.)
//		if((i.replaceControl = 1) && (std::find(controlsTriggered.begin(), controlsTriggered.end(), controlType) == controlsTriggered.end())){
//			for (auto& cTriggered:controlsTriggered){
//				std::cout<<"SM::replaceControl check. So far: "<<cTriggered<<std::endl;
//			}
//		}else{
//			std::cout<<"SM::replaceControl waitlist[controlType] output: "<<controlType<<std::endl;
//		}

		if (controlScale.compare("premises")==0){
			for (auto& f:i.premList){
				int fid = verify_premStatus(f);
				Prem_status* premToWL = changedStatus.at(fid);
				bool alreadyWaitlisted = premToWL->Prem_status::is_onWaitlist(controlType);
				if (alreadyWaitlisted == 0){
					waitlist[controlType].emplace_back(f);
					if (parameters->dangerousContacts_on==1){
						std::string fileStatus = getAny_fileStatus(premToWL);
if(verbose>1){std::cout<<"File status of new waitlist member: "<<fileStatus<<std::endl;}
						premToWL->Prem_status::set_statusWhenWaitlisted(controlType, fileStatus); // mostly interested in DC vs reported
					}
				}
			}
		} else if (controlScale.compare("premises")!=0){ // if scale is county or state
			for (auto& r:i.regionList){ // r is string of region id
				waitlistRegion[controlType].emplace_back(r);
			}
		}
	}

 	// clear newly reported
	//Start tracking numbers of newly reported premises if reports have begun.
	//Start calculating percentIncrease if number of reports have started to be tracked
	allPremReports += newPremReports.size();
	if (newPremReports.size()>0 && reportNumbers.size()==0){
		reportNumbers.emplace_back(newPremReports.size());
	}else if(newPremReports.size()>0 && reportNumbers.size()==1){ 
		reportNumbers.emplace_back(newPremReports.size());
		percentIncrease = ((reportNumbers[1]-reportNumbers[0])/reportNumbers[0])*100;
		reportNumbers.erase(reportNumbers.begin());
	}else {
		percentIncrease = newPremReports.size(); //percentIncrease is zero
		if(reportNumbers.size()>0){ //Clear the vector if there was no increase following a day of increase (starting fresh for pI calcuation)
			reportNumbers.clear();
		}
	}
	if(verbose>1){std::cout<<newPremReports.size() << " Farm(s) added. Now "<< allPremReports<< " farms have been reported."<< std::endl;}
	if(verbose>1){std::cout<<"Percent increase in number of farm reports for this time step was "<< percentIncrease<< std::endl;}


//Track the number of timesteps where there are no new prem reports or decreasing number of cases. 
	if(newPremReports.size()==0 || percentIncrease <= 0){
		reportDecreaseNumbers.emplace_back(newPremReports.size());
		decrease = reportDecreaseNumbers.size();
	}else{
		reportDecreaseNumbers.clear();
		decrease = reportDecreaseNumbers.size(); //clear vector used to track if cases uptick again.
//	std::cout<<"SM::percentIncrease calculation: clearing reportDecreaseNumbers vector"<<std::endl;
	} 
	if(verbose > 1){std::cout<<"Number of timesteps with no new reports or decreasing number of reports is: "<<decrease<<std::endl;}

	newPremReports.clear();
	newCoReports.clear();
	newStateReports.clear();

if (verbose >1){
	for (auto& c: waitlist){
	std::cout<< "On waitlist "<<c.first<<": "<<c.second.size()<<" farms"<<std::endl;
	}
	for (auto& c: waitlistRegion){
	std::cout<< "On waitlist "<<c.first<<": "<<c.second.size()<<" regions"<<std::endl;
	}
}

}
/// Passes conditions about current status of outbreak (number of premises infected, etc)
/// and checks if diagnostic rules indicate a diagnostic action should be triggered. Adds
/// premises or regions to appropriate diagnoistic waitlists.
void Status_manager::add_diagnosticWaitlistMembers(int t)
{
	
	if (verbose >1 && newPremReports.size()+newPremSuspects.size()>0){
	std::cout<<"SM::add_diagnosticWaitlistMembers: Checking for diagnostic triggers: "<<newPremReports.size()<<" new prem reports, "<<
	newPremSuspects.size()<<" new prem suspects "<<newPremVaxs.size()<<" new prem vaxs"<<std::endl;
	}
		// statsForDiagnosticTriggers in Diagnostic_manager and waitlistGroup struct defined in Control_manager

		statsForDiagnosticRules dStats{&newPremReports, &newPremSuspects, &newPremVaxs, t};
		std::vector<diagnosticWaitlistGroup> dInstructions;
		diagnosticManager->check_diagnosticRules(dStats, dInstructions); // writes waitlist info to instructions
		// add to waitlist according to instructions
		for (auto& i:dInstructions){
			std::string diagnosticType = i.diagnosticType; //stays control because it is defined that way in Control_manager
			std::string diagnosticScale = i.scaleType;
			if (diagnosticScale.compare("premises")==0){
				for (auto& f:i.premList){
					int fid = verify_premStatus(f);
					Prem_status* premToWL = changedStatus.at(fid);
					bool alreadyDiagnoisticWaitlisted = premToWL->Prem_status::is_onDiagnosticWaitlist(diagnosticType);
					if (alreadyDiagnoisticWaitlisted == 0){
						diagnosticWaitlist[diagnosticType].emplace_back(f);
						if (parameters->dangerousContacts_on==1){
								std::string fileStatus = getAny_fileStatus(premToWL);
								if(verbose>1){std::cout<<"File status of new diagnostic waitlist member: "<<fileStatus<<std::endl;}
									premToWL->Prem_status::set_statusWhenDiagnosticWaitlisted(diagnosticType, fileStatus);
								}
					}
				}
			}
		}

	// clear newly reported only if control is not on
		if (parameters->control_on==false){
			newPremReports.clear();
		}
	newPremSuspects.clear();
	newPremVaxs.clear();

if (verbose >1){
	for (auto& d: diagnosticWaitlist){
	std::cout<< "On diagnoistic waitlist "<<d.first<<": "<<d.second.size()<<" farms"<<std::endl;
	}

}

}


////For control types that have resources added during the simulation this function checks the resource levels and adds to them if necessary

/// For any control types with waitlist, checks with Control_manager to determine what
/// control is implemented for which farms/regions, based on constraint availability
/// This function also tracks how long a particular control type waitlist has remained unchanged.
void Status_manager::add_implemented(int t)
{
	for (auto& w:waitlist){ // for each control type (w.first)
		if (w.second.size()>0){
			std::vector<Farm*> updatedWaitlist;
			std::vector<Farm*> toImplement;
			controlManager->filter_constraints(w.first, w.second, updatedWaitlist, toImplement,
                                               controlResourceLevels, partiallyControlledPrems);
				
			if(w.second==updatedWaitlist){ //if the contents of the waitlist remain unchanged after filter_constraints()
				waitlistFlag[w.first]+=1;
				if(verbose>1){std::cout<<"add_implemented: no farms have come off of "<<w.first<<" waitlist this time step. Waitlist has now been halted for "<<waitlistFlag[w.first]<<"."<<std::endl;}
			} else{ //reset tracking if the waitlist changes
				waitlistFlag[w.first]=0;
				if(verbose>1){std::cout<<"add_implemented: "<<w.first<<" waitlist has changed this timestep. Setting WF to: "<<waitlistFlag[w.first]<<std::endl;}
			}

			w.second = updatedWaitlist;
			std::tuple<double, double>  duration = allControlTypes->at(w.first)->implementDuration;
 			//std::cout<<"SM:add_implemented: "<<toImplement.size()<<" is the size of to Implement at control type "<<w.first<<std::endl;

			for (auto& i:toImplement){
				// add to control list & add probPreventExposure values
				set_status(i, t, "implemented."+w.first, duration);
			}

		}
	}

	for (auto& w:waitlistRegion){ // for each control type (w.first)
		if (w.second.size()>0){
			std::vector<std::string> updatedRegionWaitlist;
			std::vector<std::string> regionsToImplement;
			controlManager->filter_constraints(w.first, w.second, updatedRegionWaitlist, regionsToImplement, &controlResourceLevels);
			w.second = updatedRegionWaitlist;
			std::tuple<double, double> duration = allControlTypes->at(w.first)->implementDuration;
			std::string regionType = allControlTypes->at(w.first)->scale;

			for (auto& i:regionsToImplement){
				// add to control list & add probPreventExposure values
				set_regionStatus(i, regionType, t, "implemented."+w.first, duration);
			}

		}
	}
}

////For diagnostic types that have resources added during the simulation this function checks the resource levels and adds to them if necessary

/// For any diagnostic types with waitlist, checks with Diagnostic_manager to determine what
/// diagnostic is started for which farms, based on constraint availability
void Status_manager::add_started(int t)
{
	toTest.clear();
	for (auto& dw:diagnosticWaitlist){ // for each diagnostic type (w.first)
		if (dw.second.size()>0){
			std::vector<Farm*> updatedDiagnosticWaitlist;
			std::vector<Farm*> toStart;


			diagnosticManager->filter_diagnosticConstraints(dw.first, dw.second, updatedDiagnosticWaitlist, toStart, diagnosticResourceLevels, partiallyTestedPrems);

			dw.second = updatedDiagnosticWaitlist;
			 std::tuple<double, double> duration = allDiagnosticTypes->at(dw.first)->startToCompleteDuration;

 			toTest[dw.first] = toStart;
 			std::cout<<"SM:add_started: "<<toTest[dw.first].size()<<" is the size of toTest at diagnostic type "<<dw.first<<std::endl;

			for (auto& i:toStart){
				// add to diagnostics list & add probPreventExposure values
				set_status(i, t, "started."+dw.first, duration);
			}
		}
	}
}

/// Returns file status of a Prem_status, or "notDangerousContact" if Prem_status does not exist
std::string Status_manager::getAny_fileStatus(Farm* f) const
{
//	std::string output;
	int fid = f->Farm::get_id();
	auto it = changedStatus.find(fid);
	if(it != changedStatus.end())
    {
        return it->second->get_fileStatus();
    }
    else
    {
        return "notDangerousContact";
    }
//	if (changedStatus.count(fid)==0){
//		output = "notDangerousContact";
//	} else {
//		output = changedStatus.at(fid)->Prem_status::get_fileStatus();
//	}
//	return output;
}

/// Returns disease status of a Prem_status, or "sus" if Prem_status does not exist
std::string Status_manager::getAny_diseaseStatus(Farm* f) const
{
	std::string output;
	int fid = f->Farm::get_id();
	if (changedStatus.count(fid)==0){
		output = "sus";
	} else {
		output = changedStatus.at(fid)->Prem_status::get_diseaseStatus();
	}
	return output;
}

/// Returns diagnostic status of a Prem_status, or "none" if Prem_status does not exist
std::string Status_manager::getAny_diagnosticStatus(Farm* f)
{
	std::string output;
	int fid = f->Farm::get_id();
	if (changedStatus.count(fid)==0){
		output = "none";
	} else {
		output = changedStatus.at(fid)->Prem_status::get_diagnosticStatuses();
	}
	//std::cout<<"SM: get_diagnosticStatus output "<<output<<std::endl;
	return output;
}

void Status_manager::get_premsWithStatus(std::vector<std::string> status_vector, std::vector<Farm*>& output)
{
    std::vector<Farm*> temp_output;
    for(std::string status : status_vector)
    {
        std::vector<Farm*> partial_output;
        get_premsWithStatus(status, partial_output);
        temp_output.insert(temp_output.end(), partial_output.begin(), partial_output.end());
    }
    output.swap(temp_output);
}

/// Copies all premises with current disease status status (as of last call to updates) into input vector
void Status_manager::get_premsWithStatus(std::string status, std::vector<Farm*>& output)
{
	if (diseaseStatuses.count(status)==1){
		std::vector<Prem_status*> outputPS; // return farms in [lo, end)
		get_premStatusesWithStatus(status, outputPS);
		std::vector<Farm*> tempOutput;
		for (auto& f:outputPS){
			tempOutput.emplace_back(allPrems->at(f->Farm::get_id()));
		}
		output.swap(tempOutput);
	}
}

void Status_manager::get_premStatusesWithStatus(std::vector<std::string>status_vector, std::vector<Prem_status*>& output)
{
    std::vector<Prem_status*> temp_output;
    for(std::string status : status_vector)
    {
        std::vector<Prem_status*> partial_output;
        get_premStatusesWithStatus(status, partial_output);
        temp_output.insert(temp_output.end(), partial_output.begin(), partial_output.end());
    }
    output.swap(temp_output);
}

void Status_manager::get_premStatusesWithStatus(std::string status, std::vector<Prem_status*>& output)
{
    	if (diseaseStatuses.count(status)==1){
            std::vector<Prem_status*>::iterator lo_it = diseaseStatuses.at(status).units.begin() + diseaseStatuses.at(status).lo;
            std::vector<Prem_status*> temp_output(lo_it, diseaseStatuses.at(status).units.end()); // return farms in [lo, end)
            output.swap(temp_output);
    	}
}


/// Copies all premises that were status s and DC at time of waitlisting
int Status_manager::count_allDCPrems(const std::string s)
{
	int dcCount=0;
	if (controlStatuses.count(s)==1){
		std::vector<Prem_status*> outputPS(controlStatuses.at(s).units.begin(), controlStatuses.at(s).units.end());
		// separate c_type from "implemented."
		std::string c_type = s;
		c_type.replace(0, c_type.find(".")+1, "");
		for (auto& f:outputPS){
			bool wasDC = f->Prem_status::was_dcWhenWaitlisted(c_type);
			if (wasDC==1){++dcCount;}
		}
	}
 return dcCount;
}

/// Copies all premises that were status s and DDC at time of waitlisting
int Status_manager::count_allDDCPrems(const std::string s)
{
	int dcCount=0;
	if (diagnosticStatuses.count(s)==1){
		std::vector<Prem_status*> outputPS(diagnosticStatuses.at(s).units.begin(), diagnosticStatuses.at(s).units.end());
		// separate d_type from "implemented."
		std::string d_type = s;
		d_type.replace(0, d_type.find(".")+1, "");
		for (auto& f:outputPS){
			bool wasDC = f->Prem_status::was_dcWhenDiagnosticWaitlisted(d_type);
			if (wasDC==1){++dcCount;}
		}
	}
 return dcCount;
}

/// Returns number of premises that have ever been disease status s
int Status_manager::get_totalPremsWithStatus(std::string s)
{
	int total = 0;
	if (s.compare("sus")==0){
		total = nPrems-1;
	} else if (diseaseStatuses.count(s)==1){
		total = diseaseStatuses.at(s).units.size();
	}
	return total;
}

/// Returns current number of premises with disease status s (as of last call to updates)
int Status_manager::numPremsWithStatus(std::string s)
{
	int total = 0;
	if (s.compare("sus")==0){
		total = nPrems - notSus.size();
	} else if (diseaseStatuses.count(s)==1){
		total = diseaseStatuses.at(s).units.size() - diseaseStatuses.at(s).lo; // number of farms between [lo, end)
	}
 return total;
}

/// Returns number of premises that have ever been file status s
int Status_manager::get_totalPremsWithFileStatus(std::string s)
{
	int total = 0;
	if (fileStatuses.count(s)==1){
		total = fileStatuses.at(s).units.size();
	}
	return total;
}


/// Returns current number of premises with file status s (as of last call to updates)
int Status_manager::numPremsWithFileStatus(std::string s)
{
	int total = 0;
	if (fileStatuses.count(s)==1){
		total = fileStatuses.at(s).units.size() - fileStatuses.at(s).lo; // number of farms between [lo, hi]
	}
 return total;
}

/// Returns current number of regions with file status s (as of last call to updates)
int Status_manager::get_numRegionsWithControlStatus(std::string s, std::string regionType)
{
	int total = 0;
	if (regionControlStatuses.count(s)>0){
		total = regionControlStatuses.at(s).units.size() - regionControlStatuses.at(s).lo;
	}
 return total;
}

/// \param[out] Vector of newly not-susceptible farms
/// Returns farms that became not-susceptible since this function was last called, and
/// resets placeholder to end to mark the beginning of the next iteration
void Status_manager::newNotSus(std::vector<Farm*>& output)
{
	// set start point to last endpoint marker (recentNotSus)
	std::vector<Farm*>::iterator start = std::next(notSus.begin(), recentNotSus); //advance iterator to recentNotSus places past begin()
	std::vector<Farm*>::iterator end = notSus.end();
	// return int to end of notSus
	std::vector<Farm*> recent (start,end);// return farms in [start, end]
	recent.swap(output);
	// move recentNotSus to past end of list (will be the beginning for next round)
	recentNotSus = notSus.size();
}

/// Get counties of seed premises
void Status_manager::get_seedCos(std::vector<std::string>& output)
{
	std::vector<std::string> fips;
	fips.reserve(seededFarms.size());
	for (auto& s:seededFarms){
		fips.emplace_back(s->get_fips());
	}
	fips.swap(output);
}
/// Add dangerousContacts outcomes to Prem_status (no direct access to Prem_status
/// from Grid_checker)
void Status_manager::add_potentialDC(Farm* f1, Farm* f2,
	std::vector<bool>& dcEvaluations)
{
	// f1 should already exist in changedStatus since it is infectious
	if(verbose>1){std::cout<<"SM::Storing "<<f2->Farm::get_id()<<" as possible DC of farm "<<f1->Farm::get_id()<<
	" with p|sus="<<dcEvaluations.at(static_cast<int>(DiseaseStatus::SUS))<<" and p|exp="<<
	dcEvaluations.at(static_cast<int>(DiseaseStatus::EXP))<<std::endl;}

	changedStatus.at(f1->Farm::get_id())->Prem_status::add_potentialDCInfo(f2, dcEvaluations);
	if(verbose>1){std::cout<<"SM add_potentialDC after call to premstatus_addpotentialDCINFO"<<std::endl;}

}

/// Formats results for summary output file
std::string Status_manager::formatRepSummary(int rep, int duration, double repTimeMS)
{
    std::cout<<"Entering formatRepSummary function...printing..."<<std::endl;

    //Get number of infectious premises and their respective counties.
    std::string inf_string = "inf";
    if(parameters->infectionType == InfectionType::BTB)
    {
        inf_string = "btb";
    }
	int nInf = 0;
	std::vector<Prem_status*> all_infectious;
	if (parameters->infectionType == InfectionType::FMD and diseaseStatuses.count(inf_string)>0)
    {
            all_infectious = diseaseStatuses.at(inf_string).units; //Pointer to the vector of all Prem_statuses with status given by inf_string
    }
    else if(parameters->infectionType == InfectionType::BTB)
    {
        all_infectious.assign(btb_all_infected.begin(), btb_all_infected.end());
    }
    nInf = all_infectious.size();

    std::unordered_set<County*> affected_counties;
    for(Prem_status* p : all_infectious)
    {
        affected_counties.insert(p->get_parent_county()); //Add all counties with infectious farms present to the set.
    }
    int nAffCounties = affected_counties.size();

	int nReported = 0;
	std::vector<Prem_status*> all_reports;
	if (fileStatuses.count("reported")>0)
    {
        all_reports = fileStatuses.at("reported").units; //Pointer to the vector of all Prem_statuses with status given by reported
    }
//    else if(parameters->infectionType == BTB)
//    {
//        all_reports.assign(btb_all_reported.begin(), btb_all_reported.end());
//    }
    nReported = all_reports.size();

    int nRecovered = 0;
    if(parameters->infectionType == InfectionType::BTB)
    {
        nRecovered = btb_all_recovered.size();
    }

	std::vector<int> controlCounts; // for each control type, will be implemented count, then effective count
	std::vector<int> diagnosticCounts; // for each diagnostic type, will be implemented count, then effective count

	if(parameters->control_on==1){
		for (auto& ct:(parameters->controlTypes)){
			int implemented = 0;
			int effective = 0;
			std::string controlStatusI = "implemented."+ct;
			if (controlStatuses.count(controlStatusI)>0){
if(verbose>1){std::cout<<"Total prems with status "<<controlStatusI<<" = "<<controlStatuses.at(controlStatusI).units.size()<<std::endl;}
				implemented = controlStatuses.at(controlStatusI).units.size();
			} else if (regionControlStatuses.count(controlStatusI)>0){
				implemented = regionControlStatuses.at(controlStatusI).units.size();
			}

			std::string controlStatusE = "effective."+ct;
			if (controlStatuses.count(controlStatusE)>0){
if(verbose>1){std::cout<<"Total prems with status "<<controlStatusE<<" = "<<controlStatuses.at(controlStatusE).units.size()<<std::endl;}
				effective = controlStatuses.at(controlStatusE).units.size();
			} else if (regionControlStatuses.count(controlStatusE)>0){
				effective = regionControlStatuses.at(controlStatusE).units.size();
			}

			bool showDcCount = 0;
			int dcCount = 0;
			for (auto&dcType:(parameters->dcControlTypes)){
				if (dcType == ct){
					showDcCount = 1;
					dcCount = count_allDCPrems("implemented."+ct);
				}
			}
			if (effective > implemented){ // only happens when lag from impl->eff is 0, so impl not recorded - fill in here
				implemented = effective;
				if (showDcCount==1){
					dcCount = count_allDCPrems("effective."+ct);
				}
			}

			controlCounts.emplace_back(implemented);
			controlCounts.emplace_back(effective);
			if (showDcCount==1){
				controlCounts.emplace_back(dcCount);
			}

		}
	}
	if(parameters->diagnostic_on==1){
		for (auto& dt:(parameters->diagnosticTypes)){
			int started = 0;
			int complete = 0;

			std::string diagnosticStatusS = "started."+dt;
			if (diagnosticStatuses.count(diagnosticStatusS)>0){
if(verbose>1){std::cout<<"Total prems with status "<<diagnosticStatusS<<" = "<<diagnosticStatuses.at(diagnosticStatusS).units.size()<<std::endl;}
				started = diagnosticStatuses.at(diagnosticStatusS).units.size();
			}

			std::string diagnosticStatusC = "complete."+dt;
			if (diagnosticStatuses.count(diagnosticStatusC)>0){
if(verbose>1){std::cout<<"Total prems with status "<<diagnosticStatusC<<" = "<<diagnosticStatuses.at(diagnosticStatusC).units.size()<<std::endl;}
				complete = diagnosticStatuses.at(diagnosticStatusC).units.size();
			}

			bool showDcCount = 0;
			int dcCount = 0;
			for (auto&dcType:(parameters->dcDiagnosticTypes)){
				if (dcType == dt){
					showDcCount = 1;
					dcCount = count_allDDCPrems("started."+dt);
				}
			}
			if (complete > started){ // only happens when lag from start->comp is 0, so start not recorded - fill in here
				started = complete;
				if (showDcCount==1){
					dcCount = count_allDDCPrems("complete."+dt);
				}
			}

			diagnosticCounts.emplace_back(started);
			diagnosticCounts.emplace_back(complete);
			if (showDcCount==1){
				diagnosticCounts.emplace_back(dcCount);
			}

		}
	}

	std::vector<int> seedIDs; seedIDs.reserve(seededFarms.size());
	for (auto& sf:seededFarms){
		seedIDs.emplace_back(sf->Farm::get_id());
	}
	std::string seeds = vecToCommaSepString(seedIDs);

	std::vector<std::string> seedFips;
	get_seedCos(seedFips);
	std::string seedCos = vecToCommaSepString(seedFips);

	double repTimeS = repTimeMS/1000;
	std::string toPrint;
	addItemTab(toPrint, rep); // rep #
if(verbose>1){std::cout<<"rep "<<rep;}
	addItemTab(toPrint, nInf); // # total infectious (includes seeds)
if(verbose>1){std::cout<<", nInf "<<nInf;}
	addItemTab(toPrint, nAffCounties); //Number of counties with infectious or exposed farms present.
if(verbose>1){std::cout<<", nAffCounties "<<nAffCounties<<std::endl;}
	addItemTab(toPrint, duration); // duration of epidemic
if(verbose>1){std::cout<<", duration "<<duration;}
	addItemTab(toPrint, seeds); // seed farm(s)
if(verbose>1){std::cout<<", seeds "<<seeds;}
	addItemTab(toPrint, seedCos); // seed county(s)
if(verbose>1){std::cout<<", seedCos "<<seedCos;}
	addItemTab(toPrint, repTimeS); // runtime
if(verbose>1){std::cout<<", repTimeSec "<<repTimeS<<std::endl;}
	addItemTab(toPrint, nReported); // # total reported premises
if(verbose>1){std::cout<<", nReported "<<nReported;}
	if (controlCounts.size()>0){
		for (auto& cc:controlCounts){
			addItemTab(toPrint, cc); // # total prems with effective control
		}
	}
	if (diagnosticCounts.size()>0){
		for (auto& dd:diagnosticCounts){
			addItemTab(toPrint, dd); // # total prems with effective diagnostics
		}
	}
	if (parameters->dangerousContacts_on==1){
		double sum = std::accumulate(dcsPerIP.begin(), dcsPerIP.end(), 0.0);
		double meanDCsPerIP = sum / dcsPerIP.size();
		if (dcsPerIP.size() == 0){
			meanDCsPerIP = 0;
		}
		addItemTab(toPrint, meanDCsPerIP);
	}

    if(parameters->infectionType == InfectionType::BTB)
    {
        addItemTab(toPrint, btb_n_freshly_infected_animals);
        addItemTab(toPrint, nRecovered);
    }

	toPrint.replace(toPrint.end()-1, toPrint.end(), "\n"); // add line break at end
if(verbose>1){std::cout<<"toPrint: "<<toPrint<<std::endl;}

	return toPrint;
}

/// Format controlSummary file. Pulls and prints each farm's
/// control status (implemented, effective, inactive) with its ID,
/// FIPS code, the timestep, and replicate. Only prints farms
/// whose controlStatuses have changed, so that not every farm
/// is printed each timestep.
std::string Status_manager::formatControlSummary(int rep, int t)
{
    std::string toPrint;
    
    if (controlStatuses.size()>0){
        for (auto& info:controlStatuses){
            std::string cStatus = std::get<0>(info); //pull control status
            std::vector<Prem_status*> prems = controlStatuses.at(cStatus).units; //vector of prems with each status
            
            for(auto& i:prems){
               std::string FIPS = i->Farm::get_parent_county()->get_id();
               int fid = i->Farm::get_id();
                
                if(std::find(printed.begin(), printed.end(), make_tuple(fid,cStatus)) != printed.end()){
                    continue; //skip loop if {fid,cStatus} pair found in vector. Checks to see if farm status has changed
                }else{ //if {fid,cStatus} pair not previously printed, go ahead and print new line.
                    addItemTab(toPrint, rep); // rep
                    addItemTab(toPrint, t); //timestep
                    addItemTab(toPrint, fid); //farm ID
                    addItemTab(toPrint, FIPS); //FIPS
                    addItemTab(toPrint, cStatus); //farm control status
                    toPrint.replace(toPrint.end()-1, toPrint.end(), "\n");
                    printed.emplace_back(std::make_tuple(fid,cStatus));
                }//checks to see if that fid has already been printed
            }
        }
    }
    if(parameters->printControl == 3){
        std::string out = Status_manager::formatWaitlistSummary(rep,t);
        toPrint+=out;
    }
    return toPrint;
}

/// Format waitlistSummary file. Pulls and prints each farm's
/// control status (implemented, effective, inactive) with its ID,
/// FIPS code, the timestep, and replicate. Only prints farms
/// whose controlStatuses have changed, so that not every farm
/// is printed each timestep.
std::string Status_manager::formatWaitlistSummary(int rep, int t)
{
    std::string toPrint;
            if(waitlist.size()>0){
             //   std::string waitlistedControl = waitlist.first;
                for (auto& info:waitlist){
                    std::string cStatus = info.first; //pull control status
                    std::vector<Farm*> waitlistedFarm = info.second;

                    for(auto& i:waitlistedFarm){
                       std::string FIPS = i->Farm::get_parent_county()->get_id();
                       int fid = i->Farm::get_id();
                        
                        if(std::find(printed.begin(), printed.end(), make_tuple(fid,"waitlist."+cStatus)) != printed.end()){
                            continue; //skip loop if {fid,cStatus} pair found in vector. Checks to see if farm status has changed
                        }else{ //if {fid,cStatus} pair not previously printed, go ahead and print new line.
                            addItemTab(toPrint, rep); // rep
                            addItemTab(toPrint, t); //timestep
                            addItemTab(toPrint, fid); //farm ID
                            addItemTab(toPrint, FIPS); //FIPS
                            addItemTab(toPrint, "waitlist."+cStatus); //farm control status
                            toPrint.replace(toPrint.end()-1, toPrint.end(), "\n");
                            printed.emplace_back(std::make_tuple(fid,"waitlist."+cStatus));
                        }//checks to see if that fid has already been printed
                    }
                }
            }
    return toPrint;
}

/// Formats details regarding exposed farms (excluding the initial seed) for output to
/// file. If a farm was exposed through multiple routes (i.e. local spread and shipping),
/// each of those exposure events are recorded on a line, even if they occurred at the
/// same time step (in fact, multiple exposures can only happen at the same timestep).
/// If exposure was prevented by control measures, the control types in effect at that
/// premises are listed on the same line under ControlPrevented.)
std::string Status_manager::formatDetails(int rep, int t)
// rep, ID, time, sourceID, method - not including initial seeds
{
//    std::cout<<"Entering formatDetail function...printing..."<<std::endl;
	std::string toPrint;

	if (sources.size()>0){
	for (auto& info:sources){
		int expPrem = std::get<0>(info)->Farm::get_id();
		int sourceID = std::get<1>(info)->Farm::get_id();
		std::string expFIPS = std::get<0>(info)->Farm::get_parent_county()->get_id();
		std::string sourceFIPS = std::get<1>(info)->Farm::get_parent_county()->get_id();
		int route = std::get<2>(info);
		std::string prevented = std::get<3>(info);

		addItemTab(toPrint, rep); // rep #
		addItemTab(toPrint, expPrem); // exposed prem ID
		addItemTab(toPrint, t); // time
		addItemTab(toPrint, sourceID); // source prem ID
		addItemTab(toPrint, route); // route of exposure: 0=local, 1=shipment
		addItemTab(toPrint, prevented);
		addItemTab(toPrint, expFIPS); //The FIPS code of the county of the exposed farm.
		addItemTab(toPrint, sourceFIPS); //The FIPS code of the county of the source farm.
		toPrint.replace(toPrint.end()-1, toPrint.end(), "\n"); // add line break at end
	}
	sources.clear();
	}
	return toPrint;
}

//returns prem staus pointer for corresponding farm id
Prem_status* Status_manager::get_correspondingPremStatus(int fid)
{
    Prem_status* p;
    //looks for fid in changedStatus
    // if not found, returns changedStatus.end()
    //if found, returns the pointer to prem status object
    if(changedStatus.find(fid)!=changedStatus.end())
    {
        p=changedStatus.at(fid);

    }else{
        std::cout<<"Prem staus for farm ID "<<fid<<" not found in changedStatus. Exiting..."<<std::endl;
        exit(EXIT_FAILURE);
    }


    return p;
}

// returns the norminf map
const std::unordered_map<std::string, std::vector<double>>& Status_manager::get_normInf_map() const
{
    return gridManager->Grid_manager::get_normInf_map();
}

const std::unordered_map<std::string, std::vector<double>>& Status_manager::get_normSus_map() const
{
    return gridManager->Grid_manager::get_normSus_map();
}

/// Update the within-herd spread for btb infection
/// \param[in] t Timestep
/// \param[in] quarter_idx Index of the current quarter
void Status_manager::btb_updateDisease(int t, int quarter_idx)
{
	if (diseaseStatuses.count("btb")>0){
        update(t, "btb", diseaseStatuses.at("btb")); // any 'btb' expiring -> 'rec'
        btb_updateWithinHerd(t, quarter_idx);
	}
}

void Status_manager::btb_updateWithinHerd(int t, int quarter_idx)
{
    if(diseaseStatuses.at("btb").lo != 0)
    {
        std::cout << "This should not happen." << std::endl;
        exit(EXIT_FAILURE);
    }
    if(diseaseStatuses.find("btb") != diseaseStatuses.end())
    {
        std::vector<Prem_status*>::iterator lo_it_btb = diseaseStatuses.at("btb").units.begin() + diseaseStatuses.at("btb").lo;
        std::vector<Prem_status*> prem_statuses_btb(lo_it_btb, diseaseStatuses.at("btb").units.end()); // return farms in [lo, end)
        std::vector<Prem_status*> retained_btb;
        retained_btb.reserve(prem_statuses_btb.size());
        for(Prem_status* ps : prem_statuses_btb)
        {
            bool still_infected = ps->btb_update_WH_spread(quarter_idx);
            if(still_infected)
            {
                retained_btb.push_back(ps);
            }
            else
            {
    //            changedStatus.erase(ps->get_id()); // I hope this is the only place where it needs to be removed from, else we'll have invalid pointers from now on...
                set_status(ps, t, "rec", pastEndTime);
    //            delete ps;
            }

        }
        diseaseStatuses.at("btb").units.swap(retained_btb);
    }

    if(diseaseStatuses.find("rec") != diseaseStatuses.end())
    {
        std::vector<Prem_status*>::iterator lo_it_rec = diseaseStatuses.at("rec").units.begin() + diseaseStatuses.at("rec").lo;
        std::vector<Prem_status*> prem_statuses_rec(lo_it_rec, diseaseStatuses.at("rec").units.end()); // return farms in [lo, end)
        std::vector<Prem_status*> retained_rec;
        retained_rec.reserve(prem_statuses_rec.size());
        for(Prem_status* ps : prem_statuses_rec)
        {
            bool not_infected = (ps->btb_get_N_infected() == 0);
            if(not_infected)
            {
                retained_rec.push_back(ps);
            }
        }
        diseaseStatuses.at("rec").units.swap(retained_rec);
    }
}

Prem_status* Status_manager::btb_make_trace(Prem_status* pst, int current_depth)
{
    const std::vector<Shipment*>& shipments = pst->get_received_inf_shipments();
    ++current_depth;
    //Go though the infected shipments backwards since, the last is the latest.
    Prem_status* final_origin = nullptr;
    //If prem 1 has sent to prem 2, and prem 2 has sent to prem 1 and neither are farms,
    //there will be a circular referencing that ends in an infinite loop if it's not
    //broken after a while. So unless the origin is found after 20 backtraces,
    // break out and consider the trace to be a failure.
    if(current_depth >= 20)
    {
        return final_origin;
    }
    for(auto rit = shipments.rbegin(); rit != shipments.rend(); ++rit)
    {
        Prem_status* nearest_origin = (*rit)->oPst;
        if(nearest_origin->get_prem_class()->tag != "Frm")
        {
            //It's not a farm, so source not found yet. Go deeper. I.e. to the prem one step further removed.
            return final_origin = btb_make_trace(nearest_origin, current_depth);
        }
        else
        {
            //It's a farm, so we've found the source.
            final_origin = nearest_origin;
            return final_origin;
        }

    }
    return final_origin;
}

void Status_manager::btb_updatePremSizeChangesQuarter(int t)
{
    std::vector<Prem_status*> retained_btb;
    for(auto& status_list_pair : diseaseStatuses)
    {
        std::string status = status_list_pair.first;
        statusList<Prem_status*>& status_list = status_list_pair.second;
        for(Prem_status* ps : status_list.units)
        {
            bool still_infected = ps->btb_update_size_change_quarter();
            //If this is a disease status there is a possiblility that this pemises will be "cured" by losing it's last infected animal.
            if(status == "btb")
            {
                if(still_infected)
                {
                    retained_btb.push_back(ps);
                }
                else
                {
                    set_status(ps, t, "rec", pastEndTime);
                }
            }
        }
    }
    diseaseStatuses.at("btb").units.swap(retained_btb);
}

void Status_manager::btb_evalExposure(int t)
{
    int localCount = 0;
	int shipCount = 0;
	std::vector<std::pair<Farm*, int>> toExpose;
	for (auto& e:exposureForEval){ // e = tuple of destination Farm, source Farm, route, true probability of infection, and pointer to shipment (=nullptr if not shipment spread)
		Farm* destination = std::get<0>(e);
		Farm* origin = std::get<1>(e);
		int route = std::get<2>(e);
		std::string prevented = "none"; // default value - no prevention, exposure occurs
		if (parameters->control_on == 0){

			if(route == 1) //Exposure through shipment
            {
                Shipment* s = std::get<4>(e);
                Prem_status* origin_prem_status = changedStatus.at(origin->Farm::get_id());
                s->oPst = origin_prem_status;
                //Sample number of animals in each  btb infection class in the shipment from the population
                //of the sender.
                std::vector<unsigned int> n_in_btb_classes(6,0);
                origin_prem_status->btb_subtract_from_inf_classes(s->volume, n_in_btb_classes); //Removes animals and updates prevalence in calsses accordingly.
                unsigned int n_infected = std::accumulate(n_in_btb_classes.begin()+1, n_in_btb_classes.end(), 0); //Add up all classes that are not susceptible.
                //To see if there were infected animals on the shipment there was previously a check against s->vol - n_in_btb_classes[0] > 0,
                //this works when the shipment vol can't be larger than the actual number of animals (i.e. farms and feedlots) but not when
                //the actual number of animals sent is truncated to be max the size of the prem so the sampled vol of the shipment according to USAMM (s->vol)
                //can be larger than the prem size (this is the case for markets).
                if(n_infected > 0) //Infection only happens if exposed/infectious animals are shipped. If only susc. then there will be no transmission.
                {
                    set_status(destination, t, "btb", pastEndTime); // disease status
                    set_status(destination, t, "exposed", pastEndTime);
                    Prem_status* dest_prem_status = changedStatus.at(destination->Farm::get_id());
                    dest_prem_status->btb_add_to_inf_classes(n_in_btb_classes); //Adds the shipped animals to the various classes and updates prevalence accordingly.
                    ++shipCount;
                    add_premSource(t, destination, origin, route, prevented);
                    dest_prem_status->add_infected_shipment(s);
                    s->n_infected = n_infected;
                }
                else
                {
                    s->n_infected = 0;
                }
            }
            else //Exposure through local or wildlife.
            {
                set_status(destination, t, "btb", pastEndTime); // disease status
                set_status(destination, t, "exposed", pastEndTime);
                Prem_status* dest_prem_status = changedStatus.at(destination->Farm::get_id());
                dest_prem_status->btb_expose_susceptibles(1); //One single susceptible animal becomes infected (if there is one).
                add_premSource(t, destination, origin, route, prevented);
                ++localCount;
            }
		} else if (parameters->control_on == 1){
			// determine if origin farm is controlled and transmits
            //Didn't implement vaccination code based on my assumption that we wont run any btb
			//vaccine runs (afaik there is no reasonable vaccine). If you need one, I think what can
			//be found in the FMD code (function evalExposure) is a good starting point.
			bool localTransmissionOccurs = eval_premTransmission(origin);
			bool localExposureOccurs = eval_premExposure(destination);

			if (!localTransmissionOccurs){
				if (verbose>1){std::cout<<"Transmission blocked."<<std::endl;}
				prevented = "src:"+vecToCommaSepString(changedStatus.at(origin->Farm::get_id())->Prem_status::get_controlStatuses()); // determining which of multiple control types could get complicated
				add_premSource(t, destination, origin, route, prevented);
			}

			if (!localExposureOccurs){
				if (verbose>1){std::cout<<"Exposure blocked."<<std::endl;}
				prevented = "exp:"+vecToCommaSepString(changedStatus.at(destination->Farm::get_id())->Prem_status::get_controlStatuses()); // determining which of multiple control types could get complicated
				add_premSource(t, destination, origin, route, prevented);
			}
			else if (localTransmissionOccurs && localExposureOccurs){
				if(route == 1) //Exposure through shipment
                {
                    Shipment* s = std::get<4>(e);
                    Prem_status* origin_prem_status = changedStatus.at(origin->Farm::get_id());
                    s->oPst = origin_prem_status;
                    //Sample number of animals in each  btb infection class in the shipment from the population
                    //of the sender.
                    std::vector<unsigned int> n_in_btb_classes(6,0);
                    origin_prem_status->btb_subtract_from_inf_classes(s->volume, n_in_btb_classes); //Removes animals and updates prevalence in calsses accordingly.
                    if(s->volume - n_in_btb_classes[0] > 0) //Infection only happens if exposed/infectious animals are shipped. If only susc. then there will be no transmission.
                    {
                        set_status(destination, t, "btb", pastEndTime); // disease status
                        set_status(destination, t, "exposed", pastEndTime);
                        Prem_status* dest_prem_status = changedStatus.at(destination->Farm::get_id());
                        dest_prem_status->btb_add_to_inf_classes(n_in_btb_classes); //Addes the shipped animals to the various classes and updates prevalence accordingly.
                        add_premSource(t, destination, origin, route, prevented);
                        ++shipCount;
                        dest_prem_status->add_infected_shipment(s);
                        int n_infected = std::accumulate(n_in_btb_classes.begin()+1, n_in_btb_classes.end(), 0);
                        s->n_infected = n_infected;
                    }
                    else
                    {
                        s->n_infected = 0;
                    }
                }
                else
                {
                    set_status(destination, t, "btb", pastEndTime); // disease status
                    set_status(destination, t, "exposed", pastEndTime);
                    Prem_status* dest_prem_status = changedStatus.at(destination->Farm::get_id());
                    dest_prem_status->btb_expose_susceptibles(1); //One single susceptible animal becomes infected (if there is one).
                    add_premSource(t, destination, origin, route, prevented);
                    ++localCount;
                }


            } // end if transmission occurs
		} // end "if control on"
	} // end "for each exposure to evaluate"

	if(verbose>0){
		std::cout<<localCount<<" local infections (may include double-counted)"<<std::endl;
		std::cout<<shipCount<<" shipping infections (may include local)"<<std::endl;
	}
	exposureForEval.clear();
}

void Status_manager::btb_updateRecords()
{
    std::vector<Prem_status*> prem_statuses; prem_statuses.reserve(100000);
    //Record new recovery events.
    get_premStatusesWithStatus("rec", prem_statuses);
    for(Prem_status* pst : prem_statuses)
    {
        btb_all_recovered.insert(pst);
    }

    prem_statuses.clear();
    get_premStatusesWithStatus("reported", prem_statuses);
    for(Prem_status* pst : prem_statuses)
    {
        btb_all_reported.insert(pst);
    }

    prem_statuses.clear();
    //Record new infections and currently infected.
    get_premStatusesWithStatus("btb", prem_statuses);
    for(Prem_status* pst : prem_statuses)
    {
        btb_n_freshly_infected_animals += pst->btb_get_new_external_infections_since_last_t() +
                                          pst->btb_get_new_WH_infections_since_last_t();
        btb_all_infected.insert(pst);
    }
}
