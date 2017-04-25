// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "OSC_Sim.h"

OSC_Sim::OSC_Sim(const Parameters_OSC_Sim& params,const int id){
    // Set parameters of Simulation base class
    initializeSimulation(params,id);
    // Set additional parameters
    // Object Parameters
    // Excitons
    Exciton_generation_rate = params.Exciton_generation_rate;
    Exciton_lifetime = params.Exciton_lifetime;
    R_exciton_hopping = params.R_exciton_hopping;
    FRET_cutoff = params.FRET_cutoff;
    // Test Parameters
    Enable_exciton_diffusion_test = params.Enable_exciton_diffusion_test;
    N_tests = params.N_tests;
    // Lattice Parameters
    Enable_gaussian_dos = params.Enable_gaussian_dos;
    Site_energy_stdev = params.Site_energy_stdev;
    Enable_exponential_dos = params.Enable_exponential_dos;
    Site_energy_urbach = params.Site_energy_urbach;
    // Output files
    //
    // Initialize energetic disorder
    if(Enable_gaussian_dos){
        site_energies.assign(getNumSites(),0);
        createGaussianDOSVector(site_energies,0,Site_energy_stdev,getId());
    }
    else if(Enable_exponential_dos){
        site_energies.assign(getNumSites(),0);
        createExponentialDOSVector(site_energies,0,Site_energy_urbach,getId());
    }
    else{
        site_energies.push_back(0);
    }
    // Initialize counters
    N_excitons = 0;
    N_excitons_created = 0;
    N_excitons_recombined = 0;
    // Initialize sites
    Site_OSC site;
    site.clearOccupancy();
    sites.assign(getNumSites(),site);
    unique_ptr<Site> site_ptr;
    for(int i=0;i<getNumSites();i++){
        if(Enable_gaussian_dos || Enable_exponential_dos){
            sites[i].setEnergyIt(site_energies.begin()+i);
        }
        else{
            sites[i].setEnergyIt(site_energies.begin());
        }
        site_ptr = unique_ptr<Site>(&sites[i]);
        addSite(site_ptr);
    }
    // Initialize exciton creation event
    Coords dest_coords;
    double R_exciton_generation = Exciton_generation_rate*getNumSites()*intpow(1e-7*getUnitSize(),3);
    exciton_creation_event.calculateEvent(dest_coords,0,0,getTemperature(),R_exciton_generation);
    unique_ptr<Event> event_ptr = unique_ptr<Event>(&exciton_creation_event);
    exciton_creation_it = addEvent(event_ptr);
}

double OSC_Sim::calculateDiffusionLength_avg(){
    return vector_avg(diffusion_distances);
}

double OSC_Sim::calculateDiffusionLength_stdev(){
    return vector_stdev(diffusion_distances);
}

Coords OSC_Sim::calculateExcitonCreationCoords(){
    Coords dest_coords;
    bool success = false;
    while(!success){
        dest_coords = getRandomCoords();
        if(loggingEnabled()){
            ostringstream msg;
            msg << "Attempting to create exciton at " << dest_coords.x << "," << dest_coords.y << "," << dest_coords.z << "." << endl;
            logMSG(msg);
        }
        if(isOccupied(dest_coords)){
            continue;
        }
        else{
            success = true;
            return dest_coords;
        }
    }
}

void OSC_Sim::calculateExcitonEvents(const list<unique_ptr<Object>>::iterator object_it){
    const auto exciton_it = getExcitonIt(*object_it);
    const Coords object_coords = exciton_it->getCoords();
    Coords dest_coords;
    double distance,E_delta;
    int index;
    // Exciton Hopping
    static Exciton_Hop event_hop;
    static const int range = ceil(FRET_cutoff/getUnitSize());
    static const int dim = (2*range+1);
    static vector<Exciton_Hop> hops_temp(dim*dim*dim,event_hop);
    static vector<bool> hops_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                if(i==0 && j==0 && k==0){
                    hops_valid[index] = false;
                    continue;
                }
                if(!isXPeriodic() && (object_coords.x+i>=getLength() || object_coords.x+i<0)){
                    hops_valid[index] = false;
                    continue;
                }
                if(!isYPeriodic() && (object_coords.y+j>=getWidth() || object_coords.y+j<0)){
                    hops_valid[index] = false;
                    continue;
                }
                if(!isZPeriodic() && (object_coords.z+k>=getHeight() || object_coords.z+k<0)){
                    hops_valid[index] = false;
                    continue;
                }
                dest_coords.x = object_coords.x+i+calculateDX(object_coords.x,i);
                dest_coords.y = object_coords.y+j+calculateDY(object_coords.y,j);
                dest_coords.z = object_coords.z+k+calculateDZ(object_coords.z,k);
                distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                if(!((distance-0.0001)>FRET_cutoff) && !((*getSiteIt(dest_coords))->isOccupied())){
                    hops_temp[index].setObjectIt(object_it);
                    E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords));
                    hops_temp[index].calculateEvent(dest_coords,distance,E_delta,getTemperature(),0);
                    hops_valid[index] = true;
                }
                else{
                    hops_valid[index] = false;
                }
            }
        }
    }
    // Exciton Recombination
    Exciton_Recombination event_recombination;
    event_recombination.setObjectIt(object_it);
    event_recombination.calculateEvent(object_coords,0,0,0,0);
    auto recombination_list_it = exciton_recombination_events.begin();
    advance(recombination_list_it,std::distance(excitons.begin(),exciton_it));
    *recombination_list_it = event_recombination;
    // Determine the fastest available hop event
    auto hop_target_it = hops_temp.end();
    for(auto hop_it=hops_temp.begin();hop_it!=hops_temp.end();++hop_it){
        if(hops_valid[std::distance(hops_temp.begin(),hop_it)] && (hop_target_it==hops_temp.end() || hop_it->getWaitTime()<hop_target_it->getWaitTime())){
            hop_target_it = hop_it;
        }
    }
    // Compare fastest hop event with recombination event to determine fastest event for this exciton
    unique_ptr<Event> event_ptr;
    if(hop_target_it->getWaitTime() < event_recombination.getWaitTime()){
        auto hop_list_it = exciton_hop_events.begin();
        advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
        *hop_list_it = *hop_target_it;
        event_ptr = unique_ptr<Event>(&(*hop_list_it));
        setEvent(exciton_it->getEventIt(),event_ptr);
    }
    else{
        event_ptr = unique_ptr<Event>(&(*recombination_list_it));
        setEvent(exciton_it->getEventIt(),event_ptr);
    }
}

void OSC_Sim::calculateExcitonListEvents(const vector<list<unique_ptr<Object>>::iterator>& exciton_it_vec){
    for(auto it=exciton_it_vec.begin();it!=exciton_it_vec.end();++it){
        calculateExcitonEvents(*it);
    }
}

bool OSC_Sim::checkFinished(){
    if(Enable_exciton_diffusion_test){
        return (N_excitons_recombined==N_tests);
    }
    cout << "Error checking simulation finish conditions." << endl;
    return true;
}

void OSC_Sim::deleteExciton(const list<Exciton>::iterator exciton_it){
    // Locate corresponding recombine event
    auto recombination_list_it = exciton_recombination_events.begin();
    advance(recombination_list_it,std::distance(excitons.begin(),exciton_it));
    // Locate corresponding hop event
    auto hop_list_it = exciton_hop_events.begin();
    advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
    // Delete exciton
    excitons.erase(exciton_it);
    // Delete exciton recombination event
    exciton_recombination_events.erase(recombination_list_it);
    // Delete exciton hop event
    exciton_hop_events.erase(hop_list_it);
}

bool OSC_Sim::executeExcitonCreation(const list<unique_ptr<Event>>::iterator event_it){
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Determine coords for the new exciton
    const Coords coords_new = calculateExcitonCreationCoords();
    // Create the new exciton and at it to the simulation
    Exciton exciton_new(getTime(),N_excitons_created+1,coords_new);
    excitons.push_back(exciton_new);
    unique_ptr<Object> object_ptr = unique_ptr<Object>(&excitons.back());
    auto object_it = addObject(object_ptr);
    // Add an empty hop and recombine event to the corresponding lists
    Exciton_Hop hop_event;
    exciton_hop_events.push_back(hop_event);
    Exciton_Recombination recombination_event;
    exciton_recombination_events.push_back(recombination_event);
    // Update exciton counters
    N_excitons_created++;
    N_excitons++;
    // Log event
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Created exciton " << exciton_new.getTag() << " at site " << coords_new.x << "," << coords_new.y << "," << coords_new.z << "." << endl;
        logMSG(msg);
    }
    // Find all nearby excitons and calculate their events
    vector<list<unique_ptr<Object>>::iterator> neighbors = findRecalcNeighbors(coords_new);
    calculateExcitonListEvents(neighbors);
    return true;
}

bool OSC_Sim::executeExcitonHop(const list<unique_ptr<Event>>::iterator event_it){
    if(isOccupied((*event_it)->getDestCoords())){
        cout << "Exciton hop cannot be executed. Destination site is already occupied." << endl;
        return false;
    }
    else{
        // Get event info
        Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
        // Update simulation time
        incrementTime((*event_it)->getWaitTime());
        // Move the exciton in the Simulation
        moveObject((*event_it)->getObjectIt(),(*event_it)->getDestCoords());
        // Log event
        if(loggingEnabled()){
            ostringstream msg;
            msg << "Exciton " << (*((*event_it)->getObjectIt()))->getTag() << " hopped to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
            logMSG(msg);
        }
        // Find all nearby excitons and calculate their events
        vector<list<unique_ptr<Object>>::iterator> neighbors = findRecalcNeighbors(coords_initial);
        vector<list<unique_ptr<Object>>::iterator> neighbors2 = findRecalcNeighbors((*event_it)->getDestCoords());
        neighbors.insert(neighbors.end(),neighbors2.begin(),neighbors2.end());
        removeObjectItDuplicates(neighbors);
        calculateExcitonListEvents(neighbors);
        return true;
    }
}

bool OSC_Sim::executeExcitonRecombine(const list<unique_ptr<Event>>::iterator event_it){
    // Get event info
    int exciton_tag = (*((*event_it)->getObjectIt()))->getTag();
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Output diffusion distance
    if(Enable_exciton_diffusion_test){
        diffusion_distances.push_back((*((*event_it)->getObjectIt()))->calculateDisplacement());
    }
    // delete exciton and its events
    deleteExciton(getExcitonIt(*((*event_it)->getObjectIt())));
    // remove the exciton from Simulation
    removeObject((*event_it)->getObjectIt());
    // Update exciton counters
    N_excitons--;
    N_excitons_recombined++;
    // Log event
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Exciton " << exciton_tag << " recombined at site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
        logMSG(msg);
    }
    // Find all nearby excitons and calculate their events
    vector<list<unique_ptr<Object>>::iterator> neighbors = findRecalcNeighbors(coords_initial);
    calculateExcitonListEvents(neighbors);
    return true;
}

bool OSC_Sim::executeNextEvent(){
    auto event_it = chooseNextEvent();
    string event_name = (*event_it)->getName();
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Executing " << event_name << " event" << endl;
        logMSG(msg);
    }
    if(event_name.compare(Exciton_Creation::name)==0){
        return executeExcitonCreation(event_it);
    }
    else if(event_name.compare(Exciton_Hop::name)==0){
        return executeExcitonHop(event_it);
    }
    else if(event_name.compare(Exciton_Recombination::name)==0){
        return executeExcitonRecombine(event_it);
    }
    else{
        //error
        cout << "Valid event not found when calling executeNextEvent" << endl;
        return false;
    }
}

list<Exciton>::iterator OSC_Sim::getExcitonIt(const unique_ptr<Object>& object_ptr){
    for(auto exciton_it=excitons.begin();exciton_it!=excitons.end();++exciton_it){
        if(object_ptr->getTag()==exciton_it->getTag()){
            return exciton_it;
        }
    }
}

vector<double> OSC_Sim::getDiffusionData(){
    return diffusion_distances;
}

int OSC_Sim::getN_excitons_created(){
    return N_excitons_created;
}

int OSC_Sim::getN_excitons_recombined(){
    return N_excitons_recombined;
}

void OSC_Sim::outputStatus(){
    cout << getId() << ": Time = " << getTime() << " seconds.\n";
    cout << getId() << ": " << N_excitons_created << " excitons have been created and " << getN_events_executed() << " events have been executed.\n";
    cout << getId() << ": There are " << N_excitons << " excitons in the lattice.\n";
    for(auto it=excitons.begin();it!=excitons.end();++it){
        cout << getId() << ": Exciton " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
    }
    cout.flush();
}

float OSC_Sim::getSiteEnergy(const Coords& coords){
    return sites[getSiteIndex(coords)].getEnergy();
}
