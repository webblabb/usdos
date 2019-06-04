#ifndef Population_manager_h
#define Population_manager_h

#include "Status_manager.h"
#include "Farm.h"

#include <iostream>

extern int verboseLevel;


class Population_manager
{
    private:
    Status_manager* statusManagerPointer;
    
    const Parameters* p;
    std::vector<std::string> speciesOnPrems; ///< List of species on all farms provided in premises file
    //std::unordered_map<std::string,double> popPars;
    

    //set up vector of birth times
    std::unordered_map<std::string,std::vector<double>> birthTimes;
    
    int verbose; ///< Can be set to override global setting for console output

    public:
    Population_manager(Status_manager*, const Parameters*);
    ~Population_manager();
    
    void set_initialFarmSize(Farm* f); // set initial farm size for farm f
    void verify_currentSize(Farm* f); //check if there is an current size stored and if there is not, then set one

//    void addBirths(Prem_status*, std::string sp);// adds animals to a farm via a birth pulse
//    void addAnimals(Farm* f, int t); // adds animals via births based on current time


};




#endif // Population_manager_h
