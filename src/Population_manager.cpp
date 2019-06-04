#include "Population_manager.h"

Population_manager::Population_manager(Status_manager* in_sm, const Parameters* p_in):
    statusManagerPointer(in_sm),
    p(p_in),
    speciesOnPrems(p->species)
    //popPars(p->popParams)
{
    
    verbose = verboseLevel;

    //load in any necessary parameters
    


    //which months are there births for beef
    std::vector<double> tvec;
    tvec.push_back(3);
    tvec.push_back(4);
    tvec.push_back(5);
    
    //fill in birth times
    for (auto &sp:speciesOnPrems){
        birthTimes[sp] = tvec;
    }
  
    
    if(verbose > 0){
        std::cout << "Population manager initiated."<<std::endl;
    }

}

Population_manager::~Population_manager()
{
}



// set initial farm size for farm f
void Population_manager::set_initialFarmSize(Farm* f){
    //get farm id
    int fid = f->Farm::get_id();
    Prem_status* prem_status_fid=statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object
    
    //for the species on the farm
    for (auto& sp:speciesOnPrems){
        double number = double(f->get_size(sp)); //find the size of each species
        //set the current size
        prem_status_fid->Prem_status::set_currentSize(sp, number);
    }
}

//check if there is an current size stored and if there is not, then set one
void Population_manager::verify_currentSize(Farm* f){
    int fid = f->Farm::get_id();
    Prem_status* prem_status_fid=statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object

    int total=0;
    for (auto& sp:speciesOnPrems){
        total +=prem_status_fid->Prem_status::get_currentSize(sp);
    }
    // this might be an issue if farms ever get to zero
    if(total==0){
        set_initialFarmSize(f);
    }
}
//
//// add animals to a farm via a birth pulse
//void Population_manager::addBirths(Prem_status* prem_status_fid, std::string sp){
//    
//        double birth=popPars[sp]; //species specific birth rate
//    
//        double N=prem_status_fid->Prem_status::get_currentSize(sp); //get the current size for species sp
//    
//        double additions=birth*N; //how many births to add this month
//    
//        double number=N+additions; //add births to the total population size
//        
//        if(verbose > 1){
//            std::cout << additions << sp << "animals added to farm"<<std::endl;
//        }
//    
//        prem_status_fid->Prem_status::set_currentSize(sp, number); //update current size in prem_status
//    
//}
//
//void Population_manager::addAnimals(Farm* f, int t){
//    
//    //get farm id
//    int fid = f->Farm::get_id();
//    Prem_status* prem_status_fid=statusManagerPointer->get_correspondingPremStatus(fid); // creates pointer to prem_status object
//    
//    //for the species on the farm
//    for (std::string sp:speciesOnPrems){
//      
//        //find the species specific time vector
//        std::vector<double> tvec=birthTimes[sp];
//        
//            if(std::find(tvec.begin(), tvec.end(), t) != tvec.end()) {
//                // tvec contains current time t
//                addBirths(prem_status_fid,sp); // add births the farm
//            }
//        }
//}

