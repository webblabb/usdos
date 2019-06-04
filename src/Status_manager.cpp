#include "Status_manager.h"

/// Establishes sequences of statuses for disease, file status, and control statuses.
/// Sets seed farm(s) as exposed, using index reporting time for reporting.
Status_manager::Status_manager(std::vector<Farm*>& focalFarms, const Parameters* parameters,
	Grid_manager* grid, Control_manager* control) :
		seededFarms(focalFarms), // saved for output
		parameters(parameters),
		controlManager(control),
		gridManager(grid),
    allPrems(grid->get_allFarms()),
    allCounties(grid->get_allCounties()),
    allStates(grid->get_allStates()),
    allControlTypes(control->get_controlTypes()),
    controlResources(control->get_controlResources()),
    controlReleaseSchedule(control->get_resourceBoostSchedule()),
    recentNotSus(0),
    pastEndTime(std::make_tuple(parameters->timesteps+100, 0)),
    nPrems(allPrems->size()),
    species(parameters->species) // store species for formatting later
{
	verbose = verboseLevel;

	// Specify duration of each disease status and what follows
	statusShift exp {parameters->latencyParams, "inf"};
	statusSequences["exp"] = exp;

	statusShift inf {parameters->infectiousParams, "imm"};
	statusSequences["inf"] = inf;

	statusShift imm {pastEndTime, "NA"};
	statusSequences["imm"] = imm;

	// set seed farms as exposed
	for (auto& f:focalFarms){
		set_status(f, 1, "exp", parameters->latencyParams); // "exp" = disease status exposed
        //get farm id
        int fid = f->Farm::get_id();
        //set exposure time for source farms
        //record exposedBy by itself
        //route is local spread
        //not blocked
        changedStatus.at(fid)->Prem_status::add_exposureSource(1, f, 0, "none");
	}

	if (parameters->control_on == true){
		statusShift exposed {pastEndTime, "reported"};
		// End times for fileStatus 'exposed' differs depending on dangerousContact status, handled in expose() and set_status(). pastEndTime put in as placeholder.
		statusSequences["exposed"] = exposed;
		statusShift reported {pastEndTime, "NA"}; // 'reported' status is permanent
		statusSequences["reported"] = reported;

		// Specify duration of each control status and what follows
		for (auto& ct:(parameters->controlTypes)){
			std::string c_type = ct;

			statusShift implemented {parameters->implementToEffectiveLag.at(c_type), "effective."+c_type};
			statusSequences["implemented."+c_type] = implemented;

			statusShift effective {parameters->effectiveToInactiveLag.at(c_type), "inactive."+c_type};
			statusSequences["effective."+c_type] = effective;

			statusShift inactive {pastEndTime, "NA"};
			statusSequences["inactive."+c_type] = inactive;
		}

		// set seed farms as exposed (fileStatus), for reporting
		std::tuple<double, double> iReport = parameters->indexReportLag;
		for (auto& f:focalFarms){
			int indexReportTime = normDelay(iReport); // current time + index report delay
			std::tuple<double, double> indexReportTuple = std::make_tuple(indexReportTime, 0);
			set_status(f, 1, "exposed", indexReportTuple); // file status
		}
	}

if (verbose>1){
	std::cout<<focalFarms.size()<<" exposed premises initiated. End times for exposed premises: "<<std::endl;
	for (auto& e : diseaseStatuses.at("exp").units){
		std::cout << e->get_end("exp") << ", ";
	}
	std::cout<<std::endl;
}

if(verbose > 0){
	std::cout << "Status manager initiated."<<std::endl;
	}
}

Status_manager::~Status_manager()
{
	for (auto& f:changedStatus){delete f.second;}
	for (auto& c:changedCoStatus){delete c.second;}
	for (auto& c:changedStateStatus){delete c.second;}
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
			status.compare("imm")==0){

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

	// if status is a file status (applies to all types of control):
	// other fileStatus: notDangerousContact
	} else if (status.compare("dangerousContact")==0 ||
						 status.compare("exposed")==0 ||
						 status.compare("reported")==0){

		if(verbose>1){std::cout<<status<<" (file status)"<<std::endl;}

		changedStatus.at(fid)->Prem_status::set_fileStatus(status);
		// add to appropriate file-status list
		if (fileStatuses.count(status)==0){
			fileStatuses[status].lo = 0; // set placeholder at beginning
		}
		fileStatuses.at(status).units.emplace_back(changedStatus.at(fid));

		if (status.compare("reported")==0){
			newPremReports.emplace_back(changedStatus.at(fid));
			if (parameters->dangerousContacts_on==1){
				// determine which of potential dangerous contacts will be dangerous contacts
				std::unordered_map<Farm*, std::unordered_map<std::string, bool>> pDCs =
					*(changedStatus.at(fid)->Prem_status::get_potentialDCs());
if(verbose>1 && pDCs.size()>0){std::cout<<"SM::Reported farm has "<<pDCs.size()<<" potential DCs"<<std::endl;}
				for (auto& dc:pDCs){ // dc.first is prem to evaluated as DC

if (verbose>1){std::cout<<dc.first->Farm::get_id()<<" is DC/reported? "<<getAny_fileStatus(dc.first)<<std::endl;}
					if (getAny_fileStatus(dc.first).compare("dangerousContact")!=0 ||
						getAny_fileStatus(dc.first).compare("reported")!=0){ // if not already reported or already a DC

						// get current disease status
						std::string dcDiseaseStatus = getAny_diseaseStatus(dc.first);

						// determine if dc should be dangerousContact
						if (dc.second.count(dcDiseaseStatus)>0){ // dc.second is map with index status, value bool
							bool isDC = dc.second.at(dcDiseaseStatus);
	if(verbose>1){std::cout<<"SM::set_status: potential DC "<<dc.first->Farm::get_id()<<" disease status = "<<dcDiseaseStatus<<", isDC = "<<isDC<<std::endl;}
							bool sourceIsF = 1; // assume source is reported farm unless changed...
							if (dcDiseaseStatus.compare("sus")!=0){ // if not susceptible, was it exposed by a different farm?
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

/// "Removes" expires farms with a given status and if applicable, sets next status.
/// \param[in] t Current timestep, used to determine if statuses have expired and set next status
/// \param[in] inStatus statusList containing vector of farms and integer placeholder
/// \param[in] status Status that applies to inStatusList
void Status_manager::update(int t, std::string status, statusList<Prem_status*>& inStatusList)
{
	// check from iterator "lo" to end, switch expired farms to before "lo"
	if (inStatusList.lo < inStatusList.units.size()){ // not all farms have expired
		std::vector<Prem_status*>::iterator lo_it = inStatusList.units.begin();
		std::advance(lo_it, inStatusList.lo); // move iterator to lo

		for (auto it = lo_it; it < inStatusList.units.end(); it++){ // check each farm* between lo and end
				if ( (*it)->Prem_status::get_end(status) == t ){ // if validity of this status expires for this farm today
				// set to next status

				if (status.compare("exposed")==0){ // specific to file status, not disease
				// report and add to waitlist for each control type
					//If we want to add a coin test to determine whether an exposed farm will be reported, then add here
					if(get_totalPremsWithFileStatus("reported")<1){
						firstRepTime = t; //if the number of reporteds is less than 1 then record the first reported time.

						if(verbose>1){
					std::cout << "SM:: First report time recorded as"<< firstRepTime << std::endl;
						}
					}
					set_status(*it, t, "reported", pastEndTime); // in expose(), status was set to expire at report lag
					report_countyAndState(*it, t); // also report at county and state level
				} else if (statusSequences.count(status)>0){ // if this status is auto-advanced to another status
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
		} }
	}
}


/// Update file and control statuses (wrapper for update() function)
/// \param[in] t Timestep
void Status_manager::updateControl(int t)
{
	// 'sus'->'exposed' transition determined by local spread, shipping, +control
	if (fileStatuses.count("exposed")>0){
		update(t, "exposed", fileStatuses.at("exposed")); // any 'exposed' expiring -> reported+waitlisted
	}
	// 'reported' is a permanent status, no updates performed

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

///Stochastically determines if the vaccinated animals that are present on one or both
///of the premises suppress the transmission event. If the suppression is successful
///and transmission is prevented whichever of the boolean reference paramters that is
///appropriate will be set to false: if the origin herd has been vaccinated, transmission
///is prevented from that herd and transPrevented is true; if the receiveing herd is
///vaccinated, exposure is prevented and expPrevented is true; if both are vaccinated
///there is no easy way to determine which vaccinateion event prevented exposure and
///both transPrevented and expPrevented are set to true.
void Status_manager::eval_premExpPrevByVaccination(Farm* ofarm, Farm* dfarm, double trueP, int t,
                                                   bool& transPrevented, bool& expPrevented)
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
    double dfarm_sus_unvaxed = dfarm->get_sus(); //Same but susceptibility of destination farm. Will always be based on full pop regardless of partial transmission or not.
    double dfarm_sus_vaxed = dfarm_sus_unvaxed;
    auto& normInf_map = get_normInf_map();
    auto& normSus_map = get_normSus_map();
    if(parameters->partial)
    {
        ofarm_inf_unvaxed = opst->get_inf_partial_as_unvaccinated(t, parameters, gridManager->get_normInf_map());
        if(ofarm_isvaxed)
        {
            ofarm_inf_vaxed = opst->get_inf_partial_as_vaccinated(t, parameters, gridManager->get_normInf_map());
        }
        else
        {
            ofarm_inf_vaxed = ofarm_inf_unvaxed;
        }
    }
    else
    {
        ofarm_inf_unvaxed = opst->get_inf();
        ofarm_inf_vaxed = 0.0;
        for(const auto& sp_N_pair : opst->get_currentSizeUnvaccinated())
        {
            int N = sp_N_pair.second;
            if(N > 0)
            {
                ofarm_inf_vaxed += normInf_map.at(sp_N_pair.first)  * std::pow(N, parameters->infExponents.at(sp_N_pair.first));
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
                dfarm_sus_vaxed += normSus_map.at(sp_N_pair.first)  * std::pow(N, parameters->susExponents.at(sp_N_pair.first));
            }
        }
    }

    double original_rate = -std::log(-(trueP-1));
    double vaccinated_rate = original_rate * ( (ofarm_inf_vaxed * dfarm_sus_vaxed) / (ofarm_inf_unvaxed * dfarm_sus_unvaxed) );
    double p_exp_vaccinated = 1 - std::exp(-vaccinated_rate); //This is the probability that the transmission event would have occurred with vaccinated animals.
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
		if (random <= pNoTransmission){  // farm does not become exposed due to control
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
		if (random <= pNotExposed){  // farm does not become exposed due to control probability
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
void Status_manager::eval_exposure(int t)
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
                    latency_for_this_exposure = origin_prem_status->get_end("exp") - t; //Latency is inherited.
                }
//                else if(origin_prem_status->get_diseaseStatus().compare('inf') == 0) //Not sure if this should be a feature. Need to decide if bringing infectious animals onto a farm leads to instant infectiousness. /Stefan
//                {
//                    latency_for_this_exposure = 0; //Instant infectiousness.
//                }
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
                eval_premExpPrevByVaccination(origin, destination, expP, t,
                                              vaxPreventsTransmission, vaxPreventsExposure);
            }

            localTransmissionOccurs *= !vaxPreventsTransmission; //Multiplying bools is equivalent to multiplying ones and zeros (same as AND).
            localExposureOccurs *= !vaxPreventsExposure;

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
//                    else if(origin_prem_status->get_diseaseStatus().compare('inf') == 0) //Not sure if this should be a feature. Need to decide if bringing infectious animals onto a farm leads to instant infectiousness. /Stefan
//                    {
//                      latency_for_this_exposure = 0; //Instant infectiousness.
//                    }
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
		Farm* destination = allPrems->at(s->destID);

		// Only evaluate exposure and control if destination is susceptible
		if (getAny_diseaseStatus(destination).compare("sus")==0){
			Farm* origin = allPrems->at(s->origID);
			bool exposeDestination = true; // default assumption, control will turn this off
			// Exposure does NOT happen if shipping bans are effective and realized:
			if (parameters->control_on == true && allControlTypes->count("shipBan")>0){
				std::string county_fips = origin->Farm::get_parent_county()->Region::get_id();
				std::string state = origin->Farm::get_parent_state()->Region::get_id();
				std::string shipBanScale = allControlTypes->at("shipBan")->scale;
				if (shipBanScale.compare("county")==0 && changedCoStatus.count(county_fips) > 0){
					double pBan = changedCoStatus.at(county_fips)->Region_status::get_probPreventExposure();
					if (pBan > 0){ // county shipBan is effective
						double random = uniform_rand();
						if (random <= pBan){ // county shipBan is realized
							s->ban = 1;
							exposeDestination = false;
if(verbose>1){std::cout<<"SM::filter_shipments: Shipment prevented by ban"<<std::endl;}
						}
					} // end "if county shipBan is effective"
				} else if (shipBanScale.compare("state")==0 && changedStateStatus.count(state)>0){
					double pBan = changedStateStatus.at(state)->Region_status::get_probPreventExposure();
					if (pBan>0){ // state shipBan is effective
						double random = uniform_rand();
						if (random <= pBan){ // state shipBan is realized
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
				add_premForEval(destination, origin, 1, 0.0); // stores info temporarily in Prem_status
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
			if (parameters->control_on == true){
			// determine time to reporting based on dangerousContact status
				if (getAny_fileStatus(f).compare("dangerousContact") == 0){
					set_status(f,t,"exposed",parameters->dcReportLag); // altered report time for DCs
				} else if (getAny_fileStatus(f).compare("notDangerousContact") == 0){
					set_status(f,t,"exposed",parameters->nonDCReportLag);
				}
			} // end "if control on"

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
	statsForControlRules stats{&newPremReports, &newCoReports, &newStateReports, t};
	std::vector<waitlistGroup> instructions;
	controlManager->check_controlRules(stats, instructions); // writes waitlist info to instructions
	// add to waitlist according to instructions
	for (auto& i:instructions){
		std::string controlType = i.controlType;
		std::string controlScale = i.scaleType;
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


////For control types that have resources added during the simulation this function checks the resource levels and adds to them if necessary

/// For any control types with waitlist, checks with Control_manager to determine what
/// control is implemented for which farms/regions, based on constraint availability
void Status_manager::add_implemented(int t)
{
	for (auto& w:waitlist){ // for each control type (w.first)
		if (w.second.size()>0){
			std::vector<Farm*> updatedWaitlist;
			std::vector<Farm*> toImplement;
			controlManager->filter_constraints(w.first, w.second, updatedWaitlist, toImplement,
																				controlResourceLevels, partiallyControlledPrems);
			w.second = updatedWaitlist;
			std::tuple<double, double> duration = allControlTypes->at(w.first)->implementDuration;

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

/// Returns file status of a Prem_status, or "notDangerousContact" if Prem_status does not exist
std::string Status_manager::getAny_fileStatus(Farm* f) const
{
	std::string output;
	int fid = f->Farm::get_id();
	if (changedStatus.count(fid)==0){
		output = "notDangerousContact";
	} else {
		output = changedStatus.at(fid)->Prem_status::get_fileStatus();
	}
	return output;
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

/// Copies all premises with current disease status s (as of last call to updates) into input vector
void Status_manager::get_premsWithStatus(std::string s, std::vector<Farm*>& output)
{
	if (diseaseStatuses.count(s)==1){
		std::vector<Farm*> tempOutput;
		std::vector<Prem_status*>::iterator lo_it = diseaseStatuses.at(s).units.begin() + diseaseStatuses.at(s).lo;
		std::vector<Prem_status*> outputPS(lo_it, diseaseStatuses.at(s).units.end()); // return farms in [lo, end)
		for (auto& f:outputPS){
			tempOutput.emplace_back(allPrems->at(f->Farm::get_id()));
		}
		output.swap(tempOutput);
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
	std::unordered_map<std::string, bool>& dcEvaluations)
{
	// f1 should already exist in changedStatus since it is infectious
	if(verbose>1){std::cout<<"SM::Storing "<<f2->Farm::get_id()<<" as possible DC of farm "<<f1->Farm::get_id()<<
	" with p|sus="<<dcEvaluations.at("sus")<<" and p|exp="<<dcEvaluations.at("exp")<<std::endl;}

	changedStatus.at(f1->Farm::get_id())->Prem_status::add_potentialDCInfo(f2, dcEvaluations);
	if(verbose>1){std::cout<<"SM add_potentialDC after call to premstatus_addpotentialDCINFO"<<std::endl;}

}

/// Formats results for summary output file
std::string Status_manager::formatRepSummary(int rep, int duration, double repTimeMS)
{
    std::unordered_set<County*> affected_counties;
    //Get number of infectious premises and their respective counties.
	int nInf = 0;
	if (diseaseStatuses.count("inf")>0)
    {
        std::vector<Prem_status*>* all_infectious = &diseaseStatuses.at("inf").units; //Pointer to the vector of all Prem_statuses with status "inf"
        nInf = all_infectious->size();
        for(Prem_status* p : *all_infectious)
        {
            affected_counties.insert(p->get_parent_county()); //Add all counties with infectious farms present to the set.
        }
    }
    int nAffCounties = affected_counties.size();

	std::vector<int> controlCounts; // for each control type, will be implemented count, then effective count

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
	if (controlCounts.size()>0){
		for (auto& cc:controlCounts){
			addItemTab(toPrint, cc); // # total prems with effective control
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
	toPrint.replace(toPrint.end()-1, toPrint.end(), "\n"); // add line break at end
if(verbose>1){std::cout<<"toPrint: "<<toPrint<<std::endl;}

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
const std::unordered_map<std::string, double>& Status_manager::get_normInf_map() const
{
    return gridManager->Grid_manager::get_normInf_map();
}

const std::unordered_map<std::string, double>& Status_manager::get_normSus_map() const
{
    return gridManager->Grid_manager::get_normSus_map();
}
