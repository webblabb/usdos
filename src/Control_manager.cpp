#include "Control_manager.h"
#include "Control_resource.h"

Control_manager::Control_manager(const Parameters* p_in, Grid_manager* G_in)
	:
	p(p_in),
	gridManager(G_in),
	allPrems(gridManager->get_allFarms())
{
	verbose = verboseLevel;
if (p->control_on == true){

		for (auto& ct:(p->controlTypes)){
		// Make a controlType struct
		std::string c_type = ct;
		if (verboseLevel>1){std::cout<< "c_type is: " << c_type<<std::endl;}

		std::string c_scale = (p->controlScales).at(c_type);
		std::tuple<double, double> eff = (p->effectiveness).at(c_type);
		std::tuple<double, double> effDuration = (p->effectiveToInactiveLag).at(c_type);

		if (verboseLevel>1){std::cout<< "C_scale is: " << c_scale<<std::endl;}

		// only continue (and make control resource) if effectiveness and duration are > 0
		if (std::get<0>(eff) > 0 && std::get<0>(effDuration)>0){
//		if (std::get<0>(effDuration)>0){
			std::string c_constraint = (p->constraintFunctions).at(c_type);
			std::tuple<double, double> impDuration = (p->implementToEffectiveLag).at(c_type);
			allControlTypes[c_type] = new controlType{c_scale, eff, c_constraint, impDuration,
				effDuration};

				if (verboseLevel>1){std::cout<< "C_constraint is: " << c_constraint<<std::endl;}

				// Make a Control_resource object for constraints other than 'noLimit'
			if (c_constraint.compare("stateSum")==0){
				read_controlLocations_stateSum(p, c_type);  //read in landfill locations
			}  else if (c_constraint.compare("nationalLimit")==0){
				read_controlResourceBoost_nationalLimit(p, c_type); // read in vaccine bank limits
				controlResources[c_type]["all"] = new Control_resource();
			}  else if (c_constraint.compare("dailyLimit")==0){
				// get daily limits (could also set up to read from external file)
				std::vector<double> dailyLimitParams = (p->constraintFuncParams).at(c_type);
				// for this setting, the same limits are used for all resources, so only one CR is created and referred to
				controlResources[c_type]["all"] = new Control_resource();
				controlResources.at(c_type).at("all")->Control_resource::set_fixedDailyLimit(std::make_tuple(dailyLimitParams.at(0), dailyLimitParams.at(1)));
			} // else if (c_constraint.compare(0,11,"inStateDistance")==0){ // if constraint type starts with "maxDistance"
// 				std::string distanceString = c_constraint.substr(12); // extract the number to use as max distance
// 				double distance = stringToNum<double>(distanceString);
// 				// where/at what scale should this distance be stored?
			//}
			if (verboseLevel>1){std::cout<< "endo if effectiveness and duration>0 " <<std::endl;}

		} // end "if effectiveness and duration>0"
		if (verboseLevel>1){std::cout<< "endo of loop for each control type " <<std::endl;}

	} // end "for each control type"

	if (verboseLevel>0){std::cout<<"Control_manager initiated."<<std::endl;}

} // end "if control is on"

}

Control_manager::~Control_manager()
{
	for (auto&p:allControlTypes){delete p.second;}
	for (auto&r:controlResources){
		for (auto&s:r.second){
			delete s.second;
		}
	}

}

/// Expects columns: id, x, y, carcass capacity mean, carcass capacity lower, carcass capacity upper
/// Reads locations from file. Saves each location as a control resource object.
void Control_manager::read_controlLocations(const Parameters* p, std::string c_type)
{
	// find index of p->constraintFuncFileTypes == "resourceLocs"
	size_t i = 0;
	bool found = 0;
	while (i < (p->constraintFuncFileTypes).at(c_type).size() && found==0){
		if ((p->constraintFuncFileTypes).at(c_type).at(i) == "resourceLocs"){
			found = 1;
		} else {
		  i++;
		}
	}
	std::string controlLocFile = (p->constraintFuncFiles).at(c_type).at(i);

	std::ifstream f(controlLocFile);
	if(!f){std::cout << "Failed to open " << controlLocFile << " Exiting..." << std::endl; exit(EXIT_FAILURE);}
	// get daily limits (could also set up to read from external file)
	std::vector<double> dailyLimitParams = (p->constraintFuncParams).at(c_type);

	if(f.is_open())
	{
	  skipBOM(f);
    if (verbose>0){std::cout << "Control locations file open, loading locations." << std::endl;}

		while(! f.eof())
		{
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			std::vector<std::string> line_vector = split(line, ' '); // separated by space

			if(! line_vector.empty()) // if line_vector has something in it
			{
				std::string id = line_vector[0];
				//std::string state = line_vector[1];
				std::string fips = line_vector[2];
				double x = stringToNum<double>(line_vector[3]);
				double y = stringToNum<double>(line_vector[4]);
				int meanCap = round(stringToNum<double>(line_vector[5]));
				//upperCap = round(stringToNum<double>(line_vector[6]));
				//lowerCap = round(stringToNum<double>(line_vector[7]));

				resourceLocation loc{x, y, fips, ""}; // leave state empty
				// need to know what cell this CR is in
	 			controlResources[c_type][id] = new Control_resource(loc, meanCap);
	 			controlResources.at(c_type).at(id)->Control_resource::set_fixedDailyLimit(std::make_tuple(dailyLimitParams.at(0), dailyLimitParams.at(1)));
			} // end if line has something in it
		} // end while not end of file
	} // if file is open

}



/// Expects columns: id, x, y, carcass capacity mean, carcass capacity lower, carcass capacity upper
/// Reads locations from file. Saves sum of all location capacities into control resource
/// objects by state.
void Control_manager::read_controlLocations_stateSum(const Parameters* p, std::string c_type)
{
	// find index of p->constraintFuncFileTypes == "resourceLocs"
	size_t i = 0;
	bool found = 0;
	while (i < (p->constraintFuncFileTypes).at(c_type).size() && found==0){
		if ((p->constraintFuncFileTypes).at(c_type).at(i) == "resourceLocs"){
			found = 1;
		} else {
		  i++;
		}
	}
	std::string controlLocFile = (p->constraintFuncFiles).at(c_type).at(i);

	std::ifstream f(controlLocFile);
	if(!f){std::cout << "Failed to open " << controlLocFile << " Exiting..." << std::endl; exit(EXIT_FAILURE);}
	// get daily limits (could also set up to read from external file)
	std::vector<double> dailyLimitParams = (p->constraintFuncParams).at(c_type);

	if(f.is_open())
	{
	  skipBOM(f);
    if (verbose>0){std::cout << "Control locations file open, loading locations." << std::endl;}

		while(! f.eof())
		{
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			std::vector<std::string> line_vector = split(line, '\t'); // separated by space

			if(! line_vector.empty()) // if line_vector has something in it
			{
				//int id = stringToNum<int>(line_vector[0]);
				std::string state = line_vector[1];
				//std::string fips = line_vector[2];
				//double x = stringToNum<double>(line_vector[3]);
				//double y = stringToNum<double>(line_vector[4]);
				int meanCap = round(stringToNum<double>(line_vector[5]));
				//upperCap = round(stringToNum<double>(line_vector[6]));
				//lowerCap = round(stringToNum<double>(line_vector[7]));

				resourceLocation loc{0,0,"",state}; // leave x, y, and fips blank
				// store resources by state
				if (controlResources.count(c_type)<1 || controlResources.at(c_type).count(state)<1){
					controlResources[c_type][state] = new Control_resource(loc, 0); // start with totallimit 0, to be overwritten
				}
				controlResources.at(c_type).at(state)->Control_resource::set_fixedDailyLimit(std::make_tuple(dailyLimitParams.at(0), dailyLimitParams.at(1)));
				controlResources.at(c_type).at(state)->Control_resource::add_capacity(meanCap);

			} // end if line has something in it
		} // end while not end of file
	} // if file is open
}

/// Expects columns: id, x, y, carcass capacity mean, carcass capacity lower, carcass capacity upper
/// Reads locations from file. Saves sum of all location capacities into control resource
/// objects by state.
void Control_manager::read_controlResourceBoost_nationalLimit(const Parameters* p, std::string c_type)
{
    if (verbose>1){std::cout << "entering controlResourceBoost_nationalLimit" << std::endl;}

	// find index of p->constraintFuncFileTypes == "resourceBoosts"
	size_t i = 0;
	bool found = 0;
	while (i < (p->constraintFuncFileTypes).at(c_type).size() && found==0){
		if ((p->constraintFuncFileTypes).at(c_type).at(i) == "resourceBoosts"){
			found = 1;
		} else {
		  i++;
		}
	}
	std::string controlBoostFile = (p->constraintFuncFiles).at(c_type).at(i);

	std::ifstream f(controlBoostFile);
	if(!f){std::cout << "CM: Control boost file not found. Exiting..." << std::endl; exit(EXIT_FAILURE);}
	// get daily limits (could also set up to read from external file)
	std::vector<double> dailyLimitParams = (p->constraintFuncParams).at(c_type);

	if(f.is_open())
	{
	  skipBOM(f);
    if (verbose>0){std::cout << "CM: Control boost file open, loading capacity." << std::endl;}

		while(! f.eof())
		{
			std::string line;
			getline(f, line); // get line from file "f", save as "line"
			std::vector<std::string> line_vector = split(line, '\t'); // separated by space

			if(! line_vector.empty()) // if line_vector has something in it
			{
				int dayAfterFirstReport = round(stringToNum<double>(line_vector[0]));
				int doses = round(stringToNum<double>(line_vector[1]));


				resourceLocation loc{0,0,"",""}; // leave x, y, fips and state blank
				if (controlResources.count(c_type)<1){
					controlResources[c_type]["all"] = new Control_resource(loc, 0); // start with totallimit 0, to be overwritten
				    if (verbose>1){std::cout << "CM: control resource for all locations initiated" << std::endl;}
				}

				if (releaseSchedule.count(c_type)<1){
					releaseSchedule[c_type][dayAfterFirstReport] = 0; // start with totallimit 0, to be overwritten
				    if (verbose>1){std::cout << "CM: release schedule for all timepoints initiated" << std::endl;}
				}

				controlResources.at(c_type).at("all")->Control_resource::set_fixedDailyLimit(std::make_tuple(dailyLimitParams.at(0), dailyLimitParams.at(1)));
							    if (verbose>0){std::cout << "CM: control resource for all locations set with fixed daily limit" << std::endl;}

				releaseSchedule.at(c_type)[dayAfterFirstReport] = doses;
			    if (verbose>0){std::cout << "CM: releaseSchedule map filled in" << std::endl;}

			} // end if line has something in it
		} // end while not end of file
	} // if file is open

}


/// Check which types of control should be applied given conditions. Returns premises and
/// regions to add to waitlists. Called from Status_manager::add_waitlistMembers().
void Control_manager::check_controlRules(statsForControlRules& input, std::vector<waitlistGroup>& toWaitlist)
{
	std::vector<waitlistGroup> output;
	for (auto& rule:(p->controlRules)){
		std::vector<Farm*> farmOutput;
		std::vector<std::string> regionOutput;
		std::string scale = "";
		if (allControlTypes.count(rule.action)>0){
			scale = allControlTypes.at(rule.action)->scale;
			if (rule.trigger.compare("newPremReportsOverX")==0){ // if trigger=="newPremReportsOverX"
	if(verbose>1){std::cout<<"CM::check_controlRules:Applying rule "<<rule.action<<" with target "<<rule.target<<std::endl;}
				apply_rule(input.newPremReports, rule, farmOutput);
			} else if (rule.trigger.compare("newRegionReportsOverX")==0){ // if trigger=="newRegionReportsOverX"
				// get the region type based on the control type
				if (scale.compare("county")==0){
					apply_regionRule(input.newCountyReports, rule, regionOutput);
				} else if (scale.compare("state")==0){
					apply_regionRule(input.newStateReports, rule, regionOutput);
				}
			} // Can also add in options for premTotalOverX,regionTotalOverX,tOverX, but need to decide targets
		}
		waitlistGroup forThisRule {rule.action, scale, farmOutput, regionOutput};
if(verbose>1){std::cout<<"CM::check_controlRules:Adding "<<farmOutput.size()<<" premises to waitlist "<<rule.action<<std::endl;}
		output.emplace_back(forThisRule);
	}
	output.swap(toWaitlist);
}

/// Determines which farms are controlled and prioritize to add to waitlist. Called from
/// Control_manager::check_controlRules().
void Control_manager::apply_rule(std::vector<Prem_status*>* reported, const controlRule& rule,
	std::vector<Farm*>& output)
{
	std::vector<Farm*> tempOutput;
	if (reported->size() > rule.threshold){
		std::vector<Farm*> input;
		// determine which farms control applies to:
		// if applying to farms in 'reported'
		if (rule.target == 0){
			for (auto& rp:(*reported)){
				Farm* f = allPrems->at(rp->Farm::get_id()); // convert Prem_status* to Farm*
				input.emplace_back(f);
			}
			prioritize(rule.priority, input, tempOutput);

		// if applying to DCs of reported
		} else if (rule.target == -1){
			for (auto& rp:(*reported)){
				std::vector<Farm*> DCs = rp->Prem_status::get_dangerousContacts();
if(verbose>1 && DCs.size()>0){std::cout<<"Farm has "<<DCs.size()<<" dangerous contacts"<<std::endl;}
				for (auto& dc:DCs){
					input.emplace_back(dc);
				}
			}
			prioritize(rule.priority, input, tempOutput);

		// if applying to neighbors within radius
		} else if (rule.target > 0){ // if any other positive number, get farms in that radius from reported, could prioritize by closest
			bool distanceRequired = 0;
			if (rule.priority.compare("closest") == 0){
				distanceRequired = 1;
			} // if not required, some neighbors in-cell may have distance = -1

			std::multimap<double, Farm*> allTargets;
			for (auto& rp:(*reported)){
				Farm* f = allPrems->at(rp->Farm::get_id()); // convert Prem_status* to Farm*
				std::multimap<double, Farm*> neighborsOfRP;
				gridManager->get_neighborsInRadius(f, rule.target, rule.radiusSquared,
																					 distanceRequired, neighborsOfRP);
				allTargets.insert(neighborsOfRP.begin(), neighborsOfRP.end()); // sorted by distance to ANY reported premises
			}
			// remove duplicates by only adding to vector if not already present
			for (auto& n:allTargets){ // n is for neighbor
				auto inVector = std::find(input.begin(), input.end(), n.second);
				if (inVector == input.end()){ // not already in vector
					input.emplace_back(n.second);
				}
			}
			prioritize(rule.priority, input, tempOutput);
		}
	}
	tempOutput.swap(output);
}

void Control_manager::apply_regionRule(std::vector<Region_status*>* reported,
	const controlRule& rule, std::vector<std::string>& output)
{
	if (reported->size() > rule.threshold){
		if (rule.target != 0){ // regional control can only apply to self
			std::cout<<"ERROR in Control_manager::apply_rule: Regions may only have targetType = 0. Exiting...";
			exit(EXIT_FAILURE);
		}

		std::vector<std::string> tempOutput;
		for (auto& rp:(*reported)){
			std::string regionID = rp->Region::get_id();
			tempOutput.emplace_back(regionID);
		}
		output.swap(tempOutput);
	}
}

/// Puts farms in order according to priority type for addition to waitlist. Called from
/// Control_manager::apply_rule(). Note that once added to the waitlist, the order is not
/// changed - the prioritization step occurs relative to premises to be added to a
/// waitlist.
/// \param[in] priorityType One of the following controlled options: earliest
/// \param[in] input Vector of Farm*s to put in order
/// \param[in] output Vector of Farm*s for output, provided as blank
void Control_manager::prioritize(std::string priorityType, std::vector<Farm*>& input,
                                 std::vector<Farm*>& output)
{
	if (output.size()>0){
		std::cout<<"ERROR: in Control_manager::prioritize(), expecting blank vector to be provided for output. Exiting..."<<std::endl;
		exit(EXIT_FAILURE);
		}

	if (priorityType.compare("earliest")==0){ // just use contents of input
		input.swap(output);
	} // can also add options for "speciesX(largest), closest" - need to decide for largest if by species or sum
}

/// Decides which farms on waitlist come off waitlist and have control implemented,
/// depending on constraint availability.
/// \param[in] controlType Type of control for which waitlist and resources will be evaluated
/// \param[in] waitlist_in Ordered vector of Farm*s on waitlist
/// \param[in] waitlist_out Provided as blank, to be filled as updated waitlist
/// \param[in] toImplement Provided as blank, to be filled as list of Farm*s for which to implement control
/// \param[in] currentLevels Pointer to control resource levels that have been altered and are being tracked by Status_manager
void Control_manager::filter_constraints(std::string controlType, std::vector<Farm*>& waitlist_in,
	std::vector<Farm*>& waitlist_out, std::vector<Farm*>& toImplement,
	std::unordered_map<Control_resource*, int>& currentLevels,
	std::unordered_map<std::string, std::unordered_map<Farm*, int>>& partiallyControlledPrems)
{
	if (waitlist_out.size()>0 || toImplement.size()>0){
		std::cout<<"ERROR: in Control_manager::filter_constraints(), expecting blank vector to be provided for output. Exiting..."<<std::endl;
		exit(EXIT_FAILURE);
	}

	if ((allControlTypes.at(controlType)->constraintType).compare("noLimit")==0){
		waitlist_in.swap(toImplement);
	}
	else if ((allControlTypes.at(controlType)->constraintType).compare("dailyLimit")==0){
		int dailyMax = controlResources.at(controlType).at("all")->Control_resource::get_dailyLimit();
		for (auto& f:waitlist_in){
			int numAnimals = f->Farm::get_size_allSpecies();
			// check if this prem has already been partially controlled for this control type
			if (partiallyControlledPrems.count(controlType)>0){
				if (partiallyControlledPrems.at(controlType).count(f)>0){
					numAnimals = partiallyControlledPrems.at(controlType).at(f); // set to remaining number of animals
				}
			}
if (verbose>1){std::cout<<"Farm "<<f->Farm::get_id()<<" with "<<numAnimals<<" animals to control by "
<<controlType<<", capacity limited to "<<dailyMax<<", ";}

			int numToControl = numAnimals;
			// reduce against daily limit if necessary
			if (numToControl > dailyMax){
				numToControl = dailyMax;
				// remaining animals not controlled
				partiallyControlledPrems[controlType][f] = numAnimals-numToControl;
			}
if (verbose>1){std::cout<<numToControl<<" to be controlled, "<<numAnimals-numToControl<<
" saved for next time."<<std::endl;}
			if (numToControl == numAnimals){ // if all animals on premises controlled this round
					toImplement.emplace_back(f);
					partiallyControlledPrems[controlType][f] = 0; // leaves a record that control was implemented for this prem
			} else { // some but not all animals controlled
					waitlist_out.emplace_back(f); // leave on waitlist (in same order as waitlist_in), next time will use reduced # of animals
			}

		}
	}
	else if ((allControlTypes.at(controlType)->constraintType).compare("stateSum")==0){
		for (auto& f:waitlist_in){
		// determine applicable resource(s)
		std::string parentState = (f->Farm::get_parent_state())->Region::get_id();
		if (controlResources.at(controlType).count(parentState)>0){
			Control_resource* cr = controlResources.at(controlType).at(parentState);
			// if control resource is not already tracked by status manager, add with initial capacity
			if (currentLevels.count(cr)<1){
				int level = controlResources.at(controlType).at(parentState)->Control_resource::get_capacity();
				currentLevels[cr] = level;
			}
			int dailyMax = cr->Control_resource::get_dailyLimit();
			// if there is available total capacity
			if (currentLevels.at(cr) > 0){
				// get number of animals on this farm
				int numAnimals = f->Farm::get_size_allSpecies();
				// check if this prem has already been partially controlled for this control type
				if (partiallyControlledPrems.count(controlType)>0){
					if (partiallyControlledPrems.at(controlType).count(f)>0){
						numAnimals = partiallyControlledPrems.at(controlType).at(f); // set to remaining number of animals
					}
				}
if (verbose>1){std::cout<<"Farm "<<f->Farm::get_id()<<" with "<<numAnimals<<" animals to control by "
<<controlType<<", ";}
				int numToControl = numAnimals;
				// reduce against daily and total limit if necessary
				int combinedLimit = dailyMax;
				if (currentLevels.at(cr) < combinedLimit){
					combinedLimit = currentLevels.at(cr); // use whichever is smaller
				}
if (verbose>1){std::cout<<"capacity limited to "<<combinedLimit<<", ";}
				if (numToControl > combinedLimit){
					numToControl = combinedLimit;
					// remaining animals not controlled
					partiallyControlledPrems[controlType][f] = numAnimals-numToControl;
				}
if (verbose>1){std::cout<<numToControl<<" to be controlled, "<<numAnimals-numToControl<<
" saved for next time."<<std::endl;}
				// what to do at the farm level
				if (numToControl == numAnimals){ // if all animals on premises controlled this round
					toImplement.emplace_back(f);
					partiallyControlledPrems[controlType][f] = 0; // leaves a record that control was implemented for this prem
				} else { // some but not all animals controlled
					waitlist_out.emplace_back(f); // leave on waitlist (in same order as waitlist_in), next time will use reduced # of animals
				}
				currentLevels.at(cr) -= numToControl;
if (verbose>1){std::cout<<"Total capacity= "<<currentLevels.at(cr)<<std::endl;}
			} else { // no capacity, farm remains on waitlist
				waitlist_out.emplace_back(f);
			}
		} else { // no resource for this state, back on waitlist (in case another resource crops up later)
			waitlist_out.emplace_back(f);
		}
	}
}	else if ((allControlTypes.at(controlType)->constraintType).compare("nationalLimit")==0){
		for (auto& f:waitlist_in){
		// determine applicable resource(s)
		if (controlResources.at(controlType).count("all")>0){
			Control_resource* cr = controlResources.at(controlType).at("all");


			// if control resource is not already tracked by status manager, add with initial capacity
			if (currentLevels.count(cr)<1){
				int level = controlResources.at(controlType).at("all")->Control_resource::get_capacity();
				currentLevels[cr] = level;

			}
			std::vector<double> dailyLimitParams = (p->constraintFuncParams).at(controlType);
			controlResources.at(controlType).at("all")->Control_resource::set_fixedDailyLimit(std::make_tuple(dailyLimitParams.at(0), dailyLimitParams.at(1)));

			int dailyMax = cr->Control_resource::get_dailyLimit();

			// if there is available total capacity
			if (currentLevels.at(cr) > 0){
				// get number of animals on this farm
				int numAnimals = f->Farm::get_size_allSpecies();

				// check if this prem has already been partially controlled for this control type
				if (partiallyControlledPrems.count(controlType)>0){
					if (partiallyControlledPrems.at(controlType).count(f)>0){
						numAnimals = partiallyControlledPrems.at(controlType).at(f); // set to remaining number of animals
					}
				}

if (verbose>1){std::cout<<"Farm "<<f->Farm::get_id()<<" with "<<numAnimals<<" animals to control by "
<<controlType<<", ";}
				int numToControl = numAnimals;
				// reduce against daily and total limit if necessary
				int combinedLimit = dailyMax;


				if (currentLevels.at(cr) < combinedLimit){
					combinedLimit = currentLevels.at(cr); // use whichever is smaller
				}
if (verbose>1){std::cout<<"capacity limited to "<<combinedLimit<<", ";}
				if (numToControl > combinedLimit){
					numToControl = combinedLimit;
					// remaining animals not controlled
					partiallyControlledPrems[controlType][f] = numAnimals-numToControl;
				}
if (verbose>1){std::cout<<numToControl<<" to be controlled, "<<numAnimals-numToControl<<
" saved for next time."<<std::endl;}
				// what to do at the farm level
				if (numToControl == numAnimals){ // if all animals on premises controlled this round
					toImplement.emplace_back(f);
					partiallyControlledPrems[controlType][f] = 0; // leaves a record that control was implemented for this prem
				} else { // some but not all animals controlled
					waitlist_out.emplace_back(f); // leave on waitlist (in same order as waitlist_in), next time will use reduced # of animals
				}
				currentLevels.at(cr) -= numToControl;
if (verbose>1){std::cout<<"Total capacity= "<<currentLevels.at(cr)<<std::endl;}
			} else { // no capacity, farm remains on waitlist
				waitlist_out.emplace_back(f);
			}
		} else { // no resource at this time, back on waitlist (in case another resource crops up later)
			waitlist_out.emplace_back(f);
		}
	}
} else { // no change to waitlist
		waitlist_in.swap(waitlist_out);
	}

}

/// Overloaded to filter region waitlists
void Control_manager::filter_constraints(std::string controlType, std::vector<std::string>& waitlist_in,
	std::vector<std::string>& waitlist_out, std::vector<std::string>& toImplement,
	std::unordered_map<Control_resource*, int>*) // both provided as blank
{
	if (waitlist_out.size()>0 || toImplement.size()>0){
		std::cout<<"ERROR: in Control_manager::filter_constraints(), expect blank vector provided for output. Exiting..."<<std::endl;
		exit(EXIT_FAILURE);
	}

	if ((allControlTypes.at(controlType)->constraintType).compare("noLimit")==0){
		waitlist_in.swap(toImplement);
	}
	else { // no change to waitlist
		waitlist_in.swap(waitlist_out);
	}

}
