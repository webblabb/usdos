//© 2024 Colorado State University
#include "Diagnostic_manager.h"
#include "Diagnostic_resource.h"

Diagnostic_manager::Diagnostic_manager(const Parameters* p_in, Grid_manager* G_in)
	:
	p(p_in),
	gridManager(G_in),
	allPrems(gridManager->get_allFarms())
{
	verbose = verboseLevel;

if (p->diagnostic_on == true){
		for (auto& dt:(p->diagnosticTypes)){
		// Make a diagnosticType struct
		std::string d_type = dt;
		if (verboseLevel>1){std::cout<< "d_type is: " << d_type<<std::endl;}

		std::string d_scale = (p->diagnosticScales).at(d_type);
		std::tuple<double, double> sens = (p->sensitivity).at(d_type);
//		std::tuple<double, double> spec = (p->specificity).at(d_type);
		if (verboseLevel>1){std::cout<< "d_scale is: " << d_scale<<std::endl;}

		// only continue (and make diagnostic resource) if sensitivity > 0
		if (std::get<0>(sens) > 0){
			std::string d_constraint = (p->diagnosticConstraintFunctions).at(d_type);
			std::tuple<double, double> invDuration = (p->testStartToCompleteLag).at(d_type);
			allDiagnosticTypes[d_type] = new diagnosticType{d_scale, sens, d_constraint, invDuration};

			if (verboseLevel>1){std::cout<< "d_constraint is: " << d_constraint<<std::endl;}

				// Make a Diagnostic_resource object for constraints other than 'noLimit'
			if (d_constraint.compare("diagnosticDailyLimit")==0){
				// get daily limits (could also set up to read from external file)
				std::vector<double> diagnosticDailyLimitParams = (p->diagnosticConstraintFuncParams).at(d_type);
				// for this setting, the same limits are used for all resources, so only one DR is created and referred to
				diagnosticResources[d_type]["all"] = new Diagnostic_resource(0);
				diagnosticResources.at(d_type).at("all")->Diagnostic_resource::set_diagnosticFixedDailyLimit(std::make_tuple(diagnosticDailyLimitParams.at(0), diagnosticDailyLimitParams.at(1)));
			} //

		} // end "if sensitivity and specificity >0"
	} // end "for each diagnostic type"

	if (verboseLevel>0){std::cout<<"Diagnostic initiated."<<std::endl;}

} // end "if diagnostic is on"

}

Diagnostic_manager::~Diagnostic_manager()
{
	for (auto&p:allDiagnosticTypes){delete p.second;}
	for (auto&r:diagnosticResources){
		for (auto&s:r.second){
			delete s.second;
		}
	}

}


/// Check which types of diagnostic should be applied given conditions. Returns premises and
/// regions to add to diagnosticWaitlists. Called from Status_manager::add_diagnosticWaitlistMembers().
void Diagnostic_manager::check_diagnosticRules(statsForDiagnosticRules& input, std::vector<diagnosticWaitlistGroup>& toDiagnosticWaitlist)
{
	std::vector<diagnosticWaitlistGroup> output;
	for (auto& diagnosticRule:(p->diagnosticRules)){
		std::vector<Farm*> farmOutput;
		std::string scale = "";
		std::string trigger = "";
		if (allDiagnosticTypes.count(diagnosticRule.action)>0){
			scale = allDiagnosticTypes.at(diagnosticRule.action)->scale;
			if(verbose>1){std::cout<<"DM::check_diagnosticRules: scale is "<<scale<<std::endl;}
			trigger = diagnosticRule.trigger;
			if(verbose>1){std::cout<<"DM::check_diagnosticRules: trigger is "<<trigger<<std::endl;}
			if(p->testSuspects_on == true){
				if(verbose>1){std::cout<<"DM::check_diagnosticRules:testSuspects true "<<std::endl;}

			if (diagnosticRule.trigger.compare("newPremSuspectsOverX")==0){ // if trigger=="newPremSuspectsOverX"
			if(verbose>1){std::cout<<"DM::check_diagnosticRules: Applying rule "<<diagnosticRule.action<<" with target "<<diagnosticRule.target<<" trigger is newPremSuspects."<<std::endl;}
				apply_diagnosticRuleToSuspects(input.newPremSuspects, diagnosticRule, farmOutput);
			} }
					if (diagnosticRule.trigger.compare("newPremReportsOverX")==0){ // if trigger=="newPremReportsOverX"
					if(verbose>1){std::cout<<"DM::check_diagnosticRules:Applying rule "<<diagnosticRule.action<<" with target "<<diagnosticRule.target<<" trigger is newPremReports."<<std::endl;}
				apply_diagnosticRule(input.newPremReports, diagnosticRule, farmOutput);
			}
			if (diagnosticRule.trigger.compare("newPremVaxsOverX")==0){ // if trigger=="newPremVaxsOverX"
			if(verbose>1){std::cout<<"DM::check_diagnosticRules:Applying rule "<<diagnosticRule.action<<" with target "<<diagnosticRule.target<<" trigger is newPremVaxs."<<std::endl;}
		apply_diagnosticRuleToVaxs(input.newPremVaxs, diagnosticRule, farmOutput);
	}

		}
		diagnosticWaitlistGroup forThisRule {diagnosticRule.action, scale, farmOutput};
		if(verbose>1){
			std::cout<<"DM::check_diagnosticRules:Adding "<<farmOutput.size()<<" premises to diagnostic waitlist "<<diagnosticRule.action<<std::endl;
		}
		output.emplace_back(forThisRule);
	}
	output.swap(toDiagnosticWaitlist);
}

/// Determines which suspected farms are tested and prioritize to add to diagnosticWaitlist. Called from
/// Diagnostic_manager::check_diagnosticRules().
void Diagnostic_manager::apply_diagnosticRuleToSuspects(std::vector<Prem_status*>* suspected, const diagnosticRule& diagnosticRule,
	std::vector<Farm*>& output)
{
	if(verbose>1){
	std::cout<<"DM::apply_diagnosticRuleToSuspects: Initiated. Rule.threshold: "<<diagnosticRule.threshold<<" diagnosticRule.target: "<<diagnosticRule.target<<std::endl;
	std::cout<<"DM apply_diagnosticRuleToSuspects: suspected size is "<<suspected->size()<<std::endl;}

	std::vector<Farm*> tempOutput;
	if (suspected->size() > diagnosticRule.threshold){
		std::vector<Farm*> input;
		// determine which farms diagnostics applies to:
		// if applying to farms in 'suspected'
		if (p->testSuspects_on == true){
			for (auto& iv:(*suspected)){
				Farm* f = allPrems->at(iv->Farm::get_id()); // convert Prem_status* to Farm*
				input.emplace_back(f);
			}
			prioritize(diagnosticRule.priority, input, tempOutput);

		// if applying to diagnostic DCs of reported premises
		}
	}
	tempOutput.swap(output);
}


/// Determines which non-suspect farms are tested and prioritize to add to diagnosticWaitlist. Called from
/// Diagnostic_manager::check_diagnosticRules().
void Diagnostic_manager::apply_diagnosticRule(std::vector<Prem_status*>* reported, const diagnosticRule& diagnosticRule,
	std::vector<Farm*>& output)
{
	if(verbose>1){std::cout<<"DM::apply_diagnosticRules: Initiated. Rule.threshold: "<<diagnosticRule.threshold<<" diagnosticRule.target: "<<diagnosticRule.target<<std::endl;}

	std::vector<Farm*> tempOutput;
	if (reported->size() > diagnosticRule.threshold){
		std::vector<Farm*> input;
		// determine which farms diagnostics applies to:
		 if (diagnosticRule.target == -1){
				for (auto& rp:(*reported)){
						std::vector<Farm*> DCs = rp->Prem_status::get_dangerousContacts();
						if(verbose>1 && DCs.size()>0){
							std::cout<<"Farm has "<<DCs.size()<<" dangerous contacts"<<std::endl;}
										for (auto& dc:DCs){
											input.emplace_back(dc);
				}
			}
			prioritize(diagnosticRule.priority, input, tempOutput);

		// if applying to neighbors within radius of a reported premises
		} else if (diagnosticRule.target > 0){ // if any other positive number, get farms in that radius from investigated, could prioritize by closest
			bool distanceRequired = 0;
			if (diagnosticRule.priority.compare("closest") == 0){
				distanceRequired = 1;
			} // if not required, some neighbors in-cell may have distance = -1

			std::multimap<double, Farm*> allDiagnosticTargets;
			for (auto& rp:(*reported)){
				Farm* f = allPrems->at(rp->Farm::get_id()); // convert Prem_status* to Farm*
				std::multimap<double, Farm*> neighborsOfRP;
				gridManager->get_neighborsInRadius(f, diagnosticRule.target, diagnosticRule.radiusSquared,
																					 distanceRequired, neighborsOfRP);
				allDiagnosticTargets.insert(neighborsOfRP.begin(), neighborsOfRP.end()); // sorted by distance to ANY reported premises
			}
			// remove duplicates by only adding to vector if not already present
			for (auto& n:allDiagnosticTargets){ // n is for neighbor
				auto inVector = std::find(input.begin(), input.end(), n.second);
				if (inVector == input.end()){ // not already in vector
					input.emplace_back(n.second);
				}
			}
			prioritize(diagnosticRule.priority, input, tempOutput);
		}
	}
	tempOutput.swap(output);
}


/// Determines which vaccinated farms are tested and prioritize to add to diagnosticWaitlist. Called from
/// Diagnostic_manager::check_diagnosticRules().
void Diagnostic_manager::apply_diagnosticRuleToVaxs(std::vector<Prem_status*>* effectiveVax, const diagnosticRule& diagnosticRule,
	std::vector<Farm*>& output)
{
	if(verbose>1){
	std::cout<<"DM::apply_diagnosticRuleToVaxs: Initiated. Rule.threshold: "<<diagnosticRule.threshold<<" diagnosticRule.target: "<<diagnosticRule.target<<std::endl;
	std::cout<<"DM apply_diagnosticRuleToVaxs: suspected size is "<<effectiveVax->size()<<std::endl;}

	std::vector<Farm*> tempOutput;
	if (effectiveVax->size() > diagnosticRule.threshold){
		std::vector<Farm*> input;
		// determine which farms diagnostics applies to:
		// if applying to farms in 'suspected'
		if (diagnosticRule.target == 0){
			for (auto& iv:(*effectiveVax)){
				Farm* f = allPrems->at(iv->Farm::get_id()); // convert Prem_status* to Farm*
				input.emplace_back(f);
			}
			prioritize(diagnosticRule.priority, input, tempOutput);

		// if applying to diagnostic DCs of reported premises
		}
	}
	tempOutput.swap(output);
}


/// Puts farms in order according to priority type for addition to diagnostic waitlist. Called from
/// Diagnostic_manager::apply_diagnosticRule(). Note that once added to the diagnostic waitlist, the order is not
/// changed - the prioritization step occurs relative to premises to be added to a
/// diagnostic waitlist.
/// \param[in] diagnosticPriorityType One of the following controlled options: earliest
/// \param[in] input Vector of Farm*s to put in order
/// \param[in] output Vector of Farm*s for output, provided as blank
void Diagnostic_manager::prioritize(std::string diagnosticPriorityType, std::vector<Farm*>& input,
                                 std::vector<Farm*>& output)
{
	if (output.size()>0){
		std::cout<<"ERROR: in Diagnostic_manager::prioritize(), expecting blank vector to be provided for output. Exiting..."<<std::endl;
		exit(EXIT_FAILURE);
		}

	if (diagnosticPriorityType.compare("earliest")==0){ // just use contents of input
		input.swap(output);
	} // can also add options for "speciesX(largest), closest" - need to decide for largest if by species or sum
}

/// Decides which farms on waitlist come off waitlist and have diagnostics implemented,
/// depending on constraint availability.
/// \param[in] diagnosticType Type of diagnostic for which waitlist and resources will be evaluated
/// \param[in] diagnosticWaitlist_in Ordered vector of Farm*s on diagnostic waitlist
/// \param[in] diagnosticWaitlist_out Provided as blank, to be filled as updated diagnostic waitlist
/// \param[in] toStart Provided as blank, to be filled as list of Farm*s for which to start diagnostic
/// \param[in] currentDiagnosticLevels Pointer to diagnostic resource levels that have been altered and are being tracked by Status_manager
void Diagnostic_manager::filter_diagnosticConstraints(std::string diagnosticType, std::vector<Farm*>& diagnosticWaitlist_in,
	std::vector<Farm*>& diagnosticWaitlist_out, std::vector<Farm*>& toStart,
	std::unordered_map<Diagnostic_resource*, int>& currentDiagnosticLevels,
	std::unordered_map<std::string, std::unordered_map<Farm*, int>>& partiallyTestedPrems)
{
	if (diagnosticWaitlist_out.size()>0 || toStart.size()>0){
		std::cout<<"ERROR: in Diagnostic_manager::filter_constraints(), expecting blank vector to be provided for output. Exiting..."<<std::endl;
		exit(EXIT_FAILURE);
	}

	if ((allDiagnosticTypes.at(diagnosticType)->diagnosticConstraintType).compare("noLimit")==0){
		diagnosticWaitlist_in.swap(toStart);
	}
	else if ((allDiagnosticTypes.at(diagnosticType)->diagnosticConstraintType).compare("diagnosticDailyLimit")==0){
		int diagnosticDailyMax = diagnosticResources.at(diagnosticType).at("all")->Diagnostic_resource::get_diagnosticDailyLimit();
		for (auto& f:diagnosticWaitlist_in){
			int numAnimals = f->Farm::get_size_allSpecies();
			// check if this prem has already been partially tested for this diagnostic type
			if (partiallyTestedPrems.count(diagnosticType)>0){
				if (partiallyTestedPrems.at(diagnosticType).count(f)>0){
					numAnimals = partiallyTestedPrems.at(diagnosticType).at(f); // set to remaining number of animals
				}
			}
			if (verbose>1)
	{std::cout<<"Farm "<<f->Farm::get_id()<<" with "<<numAnimals<<" animals to test by "
<<diagnosticType<<", capacity limited to "<<diagnosticDailyMax<<", ";}

			int numToTest = numAnimals;
			// reduce against daily limit if necessary
			if (numToTest > diagnosticDailyMax){
				numToTest = diagnosticDailyMax;
				// remaining animals not tested
				partiallyTestedPrems[diagnosticType][f] = numAnimals-numToTest;
			}
			if (verbose>1)
	{std::cout<<numToTest<<" to be tested, "<<numAnimals-numToTest<<
" saved for next time."<<std::endl;}
			if (numToTest == numAnimals){ // if all animals on premises tested this round
					toStart.emplace_back(f);
					partiallyTestedPrems[diagnosticType][f] = 0; // leaves a record that diagnostic was implemented for this prem
			} else { // some but not all animals tested
					diagnosticWaitlist_out.emplace_back(f); // leave on diagnosticWaitlist (in same order as diagnosticWaitlist_in), next time will use reduced # of animals
			}

		}
	} //could add state sum and national limit later
	else { // no change to diagnosticWaitlist
		diagnosticWaitlist_in.swap(diagnosticWaitlist_out);
	}

}

//could later add regional constraints
