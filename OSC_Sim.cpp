// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "OSC_Sim.h"

OSC_Sim::OSC_Sim(const Parameters_OPV& params,const int id){
    // Set parameters of Simulation base class
    initializeSimulation(params,id);
    // Set Additional General Parameters
    Bias = params.Bias;
    // Morphology Parameters
    Enable_neat = params.Enable_neat; // Neat takes on donor properties
    Enable_bilayer = params.Enable_bilayer;
    Thickness_donor = params.Thickness_donor; // sites
    Thickness_acceptor = params.Thickness_acceptor; // sites
    Enable_random_blend = params.Enable_random_blend;
    Acceptor_conc = params.Acceptor_conc;
    // Test Parameters
    N_tests = params.N_tests;
    Enable_exciton_diffusion_test = params.Enable_exciton_diffusion_test;
    Enable_ToF_test = params.Enable_ToF_test;
    Polaron_type = params.Polaron_type;
    ToF_initial_polarons = params.ToF_initial_polarons;
    ToF_start_time = params.ToF_start_time;
    ToF_end_time = params.ToF_end_time;
    Enable_IQE_test = params.Enable_IQE_test;
    // Exciton Parameters
    Exciton_generation_rate_donor = params.Exciton_generation_rate_donor;
    Exciton_generation_rate_acceptor = params.Exciton_generation_rate_acceptor;
    Exciton_lifetime_donor = params.Exciton_lifetime_donor; // seconds
    Exciton_lifetime_acceptor = params.Exciton_lifetime_acceptor; // seconds
    R_exciton_hopping_donor = params.R_exciton_hopping_donor;
    R_exciton_hopping_acceptor = params.R_exciton_hopping_acceptor;
    FRET_cutoff = params.FRET_cutoff;
    E_exciton_binding_donor = params.E_exciton_binding_donor;
    E_exciton_binding_acceptor = params.E_exciton_binding_acceptor;
    R_exciton_dissociation_donor = params.R_exciton_dissociation_donor;
    R_exciton_dissociation_acceptor = params.R_exciton_dissociation_acceptor;
    Exciton_dissociation_cutoff = params.Exciton_dissociation_cutoff; // nm
    // Polaron Parameters
    Enable_phase_restriction = params.Enable_phase_restriction;
    R_polaron_hopping_donor = params.R_polaron_hopping_donor;
    R_polaron_hopping_acceptor = params.R_polaron_hopping_acceptor;
    Polaron_localization_donor = params.Polaron_localization_donor; // nm^-1
    Polaron_localization_acceptor = params.Polaron_localization_acceptor; // nm^-1
    Polaron_hopping_cutoff = params.Polaron_hopping_cutoff; // nm
    Enable_miller_abrahams = params.Enable_miller_abrahams;
    Enable_marcus = params.Enable_marcus;
    Reorganization_energy_donor = params.Reorganization_energy_donor;
    Reorganization_energy_acceptor = params.Reorganization_energy_acceptor;
        // Additional Lattice Parameters
    Enable_gaussian_dos = params.Enable_gaussian_dos;
    Energy_stdev_donor = params.Energy_stdev_donor;
    Energy_stdev_acceptor = params.Energy_stdev_acceptor;
    Enable_exponential_dos = params.Enable_exponential_dos;
    Energy_urbach_donor = params.Energy_urbach_donor;
    Energy_urbach_acceptor = params.Energy_urbach_acceptor;
    // Coulomb Calculation Parameters
    Coulomb_cutoff = params.Coulomb_cutoff; // nm
    Dielectric_donor = params.Dielectric_donor;
    Dielectric_acceptor = params.Dielectric_acceptor;
    // Output files
    //
    // Initialize counters
    N_excitons = 0;
    N_excitons_created = 0;
    N_excitons_recombined = 0;
    N_excitons_dissociated = 0;
    N_excitons_quenched = 0;
    N_electrons = 0;
    N_electrons_created = 0;
    N_electrons_recombined = 0;
    N_electrons_collected = 0;
    N_holes = 0;
    N_holes_created = 0;
    N_holes_recombined = 0;
    N_holes_collected = 0;
    N_geminate_recombinations = 0;
    N_bimolecular_recombinations = 0;
    N_electron_surface_recombinations = 0;
    N_hole_surface_recombinations = 0;
    // Initialize Sites
    Site_OSC site;
    site.clearOccupancy();
    sites.assign(getNumSites(),site);
    // Initialize Film Architecture and Site Energies
    initializeArchitecture();
    // Calculate Coulomb interactions lookup table
    double avgDielectric = (Dielectric_donor+Dielectric_acceptor)/2;
    int range = ceil((Coulomb_cutoff/getUnitSize())*(Coulomb_cutoff/getUnitSize()));
    Coulomb_table.assign(range+1,0);
    for(int i=1;i<(int)Coulomb_table.size();i++){
        Coulomb_table[i] = ((Coulomb_constant*Elementary_charge)/avgDielectric)/(1e-9*getUnitSize()*sqrt((double)i));
    }
    // Initialize exciton creation event
    if(Enable_exciton_diffusion_test || Enable_IQE_test){
        Coords dest_coords;
        R_exciton_generation_donor = Exciton_generation_rate_donor*N_donor_sites*intpow(1e-7*getUnitSize(),3);
        R_exciton_generation_acceptor = Exciton_generation_rate_acceptor*N_acceptor_sites*intpow(1e-7*getUnitSize(),3);
        exciton_creation_event.calculateEvent(dest_coords,0,0,getTemperature(),R_exciton_generation_donor+R_exciton_generation_acceptor);
        exciton_creation_it = addEvent(&exciton_creation_event);
    }
    else if(Enable_ToF_test){
        // Create initial test polarons
        generateToFPolarons();
    }
}

double OSC_Sim::calculateCoulomb(const list<unique_ptr<Object>>::iterator object_it,const Coords& coords){
    bool charge;
    if((*object_it)->getName().compare(Polaron::name)==0){
        charge = static_cast<Polaron*>(object_it->get())->getCharge();
    }
    else{
        cout << "Error! Coulomb interactions calculated for an object that is not a Polaron." << endl;
    }
    static const double avgDielectric = (Dielectric_donor+Dielectric_acceptor)/2;
    static const double image_interactions = (Elementary_charge/(16*Pi*avgDielectric*Vacuum_permittivity))*1e9;
    double Energy = 0;
    double distance;
    int distance_sq_lat;
    static const int range = ceil((Coulomb_cutoff/getUnitSize())*(Coulomb_cutoff/getUnitSize()));
    // Loop through electrons
    for(auto it=electrons.begin();it!=electrons.end();++it){
        distance_sq_lat = calculateLatticeDistanceSquared(coords,*it);
        if(!distance_sq_lat>range){
            if(!charge && it->getTag()!=(*object_it)->getTag()){
                Energy += Coulomb_table[distance_sq_lat];
            }
            else{
                Energy -= Coulomb_table[distance_sq_lat];
            }
        }
    }
    // Loop through holes
    for(auto it=holes.begin();it!=holes.end();++it){
        distance_sq_lat = calculateLatticeDistanceSquared(coords,*it);
        if(!distance_sq_lat>range){
            if(charge && it->getTag()!=(*object_it)->getTag()){
                Energy += Coulomb_table[distance_sq_lat];
            }
            else{
                Energy -= Coulomb_table[distance_sq_lat];
            }
        }
    }
    // Add electrode image charge interactions
    if(!isZPeriodic()){
        distance = getUnitSize()*((double)(getHeight()-coords.z)-0.5);
        if(!(distance-0.0001)>Coulomb_cutoff){
            Energy -= image_interactions/distance;
        }
        distance = getUnitSize()*((double)(coords.z+1)-0.5);
        if(!(distance-0.0001)>Coulomb_cutoff){
            Energy -= image_interactions/distance;
        }
    }
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
    // Exciton hopping and dissociation
    static const int range = ceil(FRET_cutoff/getUnitSize());
    static const int dim = (2*range+1);
    static vector<Exciton_Hop> hops_temp(dim*dim*dim);
    static vector<bool> hops_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                if(!checkMoveEventValidity(object_coords,i,j,k)){
                    hops_valid[index] = false;
                    continue;
                }
                dest_coords = calculateDestinationCoords(object_coords,i,j,k);
                if(isOccupied(dest_coords)){
                    hops_valid[index] = false;
                    continue;
                }
                distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                if(!((distance-0.0001)>FRET_cutoff)){
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
    auto recombination_event_it = exciton_recombination_events.begin();
    advance(recombination_event_it,std::distance(excitons.begin(),exciton_it));
    recombination_event_it->calculateEvent(object_coords,0,0,0,0);
    // Determine the fastest valid hop event
    bool No_hops_valid = true;
    auto hop_target_it = hops_temp.end();
    for(auto hop_it=hops_temp.begin();hop_it!=hops_temp.end();++hop_it){
        if(hops_valid[std::distance(hops_temp.begin(),hop_it)] && (hop_target_it==hops_temp.end() || hop_it->getWaitTime()<hop_target_it->getWaitTime())){
            hop_target_it = hop_it;
            No_hops_valid = false;
        }
    }
    // Compare fastest hop event with recombination event to determine fastest event for this exciton
    unique_ptr<Event> event_ptr;
    if(!No_hops_valid && hop_target_it->getWaitTime() < recombination_event_it->getWaitTime()){
        auto hop_list_it = exciton_hop_events.begin();
        advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
        *hop_list_it = *hop_target_it;
        setEvent(exciton_it->getEventIt(),&(*hop_list_it));
    }
    else{
        setEvent(exciton_it->getEventIt(),&(*recombination_event_it));
    }
}

void OSC_Sim::calculateObjectListEvents(const vector<list<unique_ptr<Object>>::iterator>& object_it_vec){
    for(auto it=object_it_vec.begin();it!=object_it_vec.end();++it){
        // If object is exciton
        if((**it)->getName().compare(Exciton::name)==0){
            calculateExcitonEvents(*it);
        }
        // If object is polaron
        else if((**it)->getName().compare(Polaron::name)==0){
            calculatePolaronEvents(*it);
        }
    }
}

void OSC_Sim::calculatePolaronEvents(const list<unique_ptr<Object>>::iterator object_it){
    const auto polaron_it = getPolaronIt(*object_it);
    const Coords object_coords = polaron_it->getCoords();
    Coords dest_coords;
    double distance,E_delta;
    int index;
    double Coulomb_i = calculateCoulomb(object_it,object_coords);
    // Calculate Polaron hopping and recombination events
    static const int range = ceil(Polaron_hopping_cutoff/getUnitSize());
    static const int dim = (2*range+1);
    static vector<Polaron_Hop> hops_temp(dim*dim*dim);
    static vector<Polaron_Recombination> recombinations_temp(dim*dim*dim);
    static vector<bool> events_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                if(!checkMoveEventValidity(object_coords,i,j,k)){
                    events_valid[index] = false;
                    continue;
                }
                dest_coords = calculateDestinationCoords(object_coords,i,j,k);
                // If destination site is occupied by a hole Polaron and the main Polaron is an electron, check for a possible recombination event
                if(isOccupied(dest_coords) && !polaron_it->getCharge() && siteContainsHole(dest_coords)){
                    distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                    if(!((distance-0.0001)>Polaron_hopping_cutoff)){
                        recombinations_temp[index].setObjectIt(object_it);
                        // Total energy change is sum of site energy change, Coulomb potential change, and E_potential change
                        E_delta = 0;
                        if(getSiteType(object_coords)==(short)1){
                            recombinations_temp[index].calculateEvent(dest_coords,distance,E_delta,getTemperature(),R_polaron_hopping_donor);
                        }
                        else if(getSiteType(object_coords)==(short)2){
                            recombinations_temp[index].calculateEvent(dest_coords,distance,E_delta,getTemperature(),R_polaron_hopping_acceptor);
                        }
                        events_valid[index] = true;
                    }
                }
                // If destination site is unoccupied and either phase restriction is disabled or the starting site and destination sites have the same type, check for a possible hop event
                else if(!isOccupied(dest_coords) && (!Enable_phase_restriction || getSiteType(object_coords)==getSiteType(dest_coords))){
                    distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                    if(!((distance-0.0001)>Polaron_hopping_cutoff)){
                        hops_temp[index].setObjectIt(object_it);
                        E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords))+(calculateCoulomb(object_it,dest_coords)-Coulomb_i)+(E_potential[dest_coords.z]-E_potential[object_coords.z]);
                        if(getSiteType(object_coords)==(short)1){
                            hops_temp[index].calculateEvent(dest_coords,distance,E_delta,getTemperature(),R_polaron_hopping_donor);
                        }
                        else if(getSiteType(object_coords)==(short)2){
                            hops_temp[index].calculateEvent(dest_coords,distance,E_delta,getTemperature(),R_polaron_hopping_acceptor);
                        }
                        events_valid[index] = true;
                    }
                }
                events_valid[index] = false;
            }
        }
    }
    // Calculate possible extraction event
    // Electrons are extracted at the bottom of the lattice (z=-1)
    // Holes are extracted at the top of the lattice (z=Height)
    bool No_extraction_valid = true;
    list<Polaron_Extraction>::iterator extraction_event_it;
    // If electron, charge is false
    if(!polaron_it->getCharge()){
        distance = getUnitSize()*((double)(object_coords.z+1)-0.5);
        if(!(distance-0.0001)>Polaron_hopping_cutoff){
            extraction_event_it = electron_extraction_events.begin();
            advance(extraction_event_it,std::distance(electrons.begin(),polaron_it));
            No_extraction_valid = false;
        }
    }
    // If hole, charge is true
    else{
        distance = getUnitSize()*((double)(getHeight()-object_coords.z)-0.5);
        if(!(distance-0.0001)>Polaron_hopping_cutoff){
            extraction_event_it = hole_extraction_events.begin();
            advance(extraction_event_it,std::distance(holes.begin(),polaron_it));
            No_extraction_valid = false;
        }
    }
    if(!No_extraction_valid){
        if(getSiteType(object_coords)==(short)1){
            extraction_event_it->calculateEvent(object_coords,distance,0,0,R_polaron_hopping_donor);
        }
        else if(getSiteType(object_coords)==(short)2){
            extraction_event_it->calculateEvent(object_coords,distance,0,0,R_polaron_hopping_acceptor);
        }
    }
    // Determine the fastest hop event
    bool No_hops_valid = true;
    auto hop_target_it = hops_temp.end();
    for(auto it=hops_temp.begin();it!=hops_temp.end();++it){
        if(events_valid[std::distance(hops_temp.begin(),it)] && (hop_target_it==hops_temp.end() || it->getWaitTime()<hop_target_it->getWaitTime())){
            hop_target_it = it;
            No_hops_valid = false;
        }
    }
    // Determine the fastest recombination event
    bool No_recombinations_valid = true;
    auto recombination_target_it = recombinations_temp.end();
    for(auto it=recombinations_temp.begin();it!=recombinations_temp.end();++it){
        if(events_valid[std::distance(recombinations_temp.begin(),it)] && (recombination_target_it==recombinations_temp.end() || it->getWaitTime()<recombination_target_it->getWaitTime())){
            recombination_target_it = it;
            No_recombinations_valid = false;
        }
    }
    // If no valid event is found, return
    if(No_extraction_valid && No_hops_valid && No_recombinations_valid){
        return;
    }
    // Determine the fastest valid event
    short selection = -1;
    double best_time = -1;
    if(!No_extraction_valid){
        selection = 1;
        best_time = extraction_event_it->getWaitTime();
    }
    if(!No_hops_valid && (best_time<0 || hop_target_it->getWaitTime()<best_time)){
        selection = 2;
        best_time = hop_target_it->getWaitTime();
    }
    if(!No_recombinations_valid && (best_time<0 || recombination_target_it->getWaitTime()<best_time)){
        selection = 3;
    }
    switch(selection){
        // Extraction event is fastest
        case 1:{
            setEvent(polaron_it->getEventIt(),&(*extraction_event_it));
            break;
        }
        // Hop event is fastest
        case 2:{
            list<Polaron_Hop>::iterator hop_list_it;
            // If electron, charge is false
            if(!polaron_it->getCharge()){
                hop_list_it = electron_hop_events.begin();
                advance(hop_list_it,std::distance(electrons.begin(),polaron_it));
            }
            // If hole, charge is true
            else{
                hop_list_it = hole_hop_events.begin();
                advance(hop_list_it,std::distance(holes.begin(),polaron_it));
            }
            *hop_list_it = *hop_target_it;
            setEvent(polaron_it->getEventIt(),&(*hop_list_it));
            break;
        }
        // Recombination event is fastest
        case 3:{
            list<Polaron_Recombination>::iterator recombination_list_it;
            // If electron, charge is false
            if(!polaron_it->getCharge()){
                recombination_list_it = polaron_recombination_events.begin();
                advance(recombination_list_it,std::distance(electrons.begin(),polaron_it));
            }
            // If hole, charge is true
            else{
                cout << "Error! Only electrons can initiate polaron recombination." << endl;
                return;
            }
            *recombination_list_it = *recombination_target_it;
            setEvent(polaron_it->getEventIt(),&(*recombination_list_it));
            break;
        }
        // No valid event found
        default:{
            return;
        }
    }
}

bool OSC_Sim::checkFinished(){
    if(Enable_exciton_diffusion_test){
        return (N_excitons_recombined==N_tests);
    }
    cout << "Error checking simulation finish conditions." << endl;
    return true;
}

void OSC_Sim::deleteObject(const list<unique_ptr<Object>>::iterator object_it){
    if((*object_it)->getName().compare(Exciton::name)==0){
        auto exciton_it = getExcitonIt(*object_it);
        // Locate corresponding recombination event
        auto recombination_list_it = exciton_recombination_events.begin();
        advance(recombination_list_it,std::distance(excitons.begin(),exciton_it));
        // Locate corresponding hop event
        auto hop_list_it = exciton_hop_events.begin();
        advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
        // Locate corresponding dissociation event
        auto dissociation_list_it = exciton_dissociation_events.begin();
        advance(dissociation_list_it,std::distance(excitons.begin(),exciton_it));
        // Delete exciton
        excitons.erase(exciton_it);
        // Delete exciton recombination event
        exciton_recombination_events.erase(recombination_list_it);
        // Delete exciton hop event
        exciton_hop_events.erase(hop_list_it);
        // Delete exciton dissociation event
        exciton_dissociation_events.erase(dissociation_list_it);
    }
    else if((*object_it)->getName().compare(Polaron::name)==0){
        // Electron
        if(!getObjectCharge(object_it)){

        }
        // Hole
        else{

        }
    }
    // Remove the object from Simulation
    removeObject(object_it);
}

bool OSC_Sim::executeExcitonCreation(const list<unique_ptr<Event>>::iterator event_it){
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Determine coordinates for the new exciton
    const Coords coords_new = calculateExcitonCreationCoords();
    // Create the new exciton and add it to the simulation
    Exciton exciton_new(getTime(),N_excitons_created+1,coords_new);
    excitons.push_back(exciton_new);
    auto object_it = addObject(&excitons.back());
    // Add placeholder events to the corresponding lists
    Exciton_Hop hop_event;
    exciton_hop_events.push_back(hop_event);
    Exciton_Recombination recombination_event;
    recombination_event.setObjectIt(object_it);
    exciton_recombination_events.push_back(recombination_event);
    Exciton_Dissociation dissociation_event;
    exciton_dissociation_events.push_back(dissociation_event);
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
    auto neighbors = findRecalcNeighbors(coords_new);
    calculateObjectListEvents(neighbors);
    return true;
}

bool OSC_Sim::executeExcitonDissociation(const list<unique_ptr<Event>>::iterator event_it){
    // Get event info
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    Coords coords_dest = (*event_it)->getDestCoords();
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Delete exciton and its events
    deleteObject((*event_it)->getObjectIt());
    // Generate new electron and hole
    int tag = (N_electrons_created>N_holes_created) ? (N_electrons_created+1) : (N_holes_created+1);
    if(getSiteType(coords_dest)==(short)2){
        generateHole(coords_initial,tag);
        generateElectron(coords_dest,tag);
    }
    else{
        generateElectron(coords_initial,tag);
        generateHole(coords_dest,tag);
    }
    // Update counters
    N_excitons--;
    N_excitons_dissociated++;
    // Find all nearby objects
    auto neighbors = findRecalcNeighbors(coords_initial);
    auto neighbors2 = findRecalcNeighbors(coords_dest);
    neighbors.insert(neighbors.end(),neighbors2.begin(),neighbors2.end());
    removeObjectItDuplicates(neighbors);
    // Calculate events for all nearby objects
    calculateObjectListEvents(neighbors);
    return true;
}

bool OSC_Sim::executeExcitonHop(const list<unique_ptr<Event>>::iterator event_it){
    if(isOccupied((*event_it)->getDestCoords())){
        cout << "Exciton hop cannot be executed. Destination site is already occupied." << endl;
        return false;
    }
    else{
        if(loggingEnabled()){
            ostringstream msg;
            msg << "Exciton " << (*((*event_it)->getObjectIt()))->getTag() << " hopping to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
            logMSG(msg);
        }
        return executeObjectHop(event_it);
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
    deleteObject((*event_it)->getObjectIt());
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
    auto neighbors = findRecalcNeighbors(coords_initial);
    calculateObjectListEvents(neighbors);
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
    else if(event_name.compare(Exciton_Dissociation::name)==0){
        return executeExcitonDissociation(event_it);
    }
//    else if(event_name.compare(Exciton_Intersystem_Crossing::name)==0){
//        return executeExcitonIntersystemCrossing(event_it);
//    }
//    else if(event_name.compare(Exciton_Exciton_Annihilation::name)==0){
//        return executeExcitonExcitonAnnihilation(event_it);
//    }
//    else if(event_name.compare(Exciton_Polaron_Annihilation::name)==0){
//        return executeExcitonPolaronAnnihilation(event_it);
//    }
    else if(event_name.compare(Polaron_Hop::name)==0){
        return executePolaronHop(event_it);
    }
    else if(event_name.compare(Polaron_Recombination::name)==0){
        return executePolaronRecombination(event_it);
    }
    else if(event_name.compare(Polaron_Extraction::name)==0){
        return executePolaronExtraction(event_it);
    }
    else{
        //error
        cout << "Valid event not found when calling executeNextEvent" << endl;
        return false;
    }
}

bool OSC_Sim::executeObjectHop(const list<unique_ptr<Event>>::iterator event_it){
    // Get event info
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Move the exciton in the Simulation
    moveObject((*event_it)->getObjectIt(),(*event_it)->getDestCoords());
    // Find all nearby objects
    auto neighbors = findRecalcNeighbors(coords_initial);
    auto neighbors2 = findRecalcNeighbors((*event_it)->getDestCoords());
    neighbors.insert(neighbors.end(),neighbors2.begin(),neighbors2.end());
    removeObjectItDuplicates(neighbors);
    // Calculate events for all nearby objects
    calculateObjectListEvents(neighbors);
    return true;
}

bool OSC_Sim::executePolaronExtraction(const list<unique_ptr<Event>>::iterator event_it){
    // Get event info
    bool charge = getObjectCharge((*event_it)->getObjectIt());
    int polaron_tag = (*((*event_it)->getObjectIt()))->getTag();
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Save transit time
    if(Enable_ToF_test){
        transit_times.push_back(getTime()-(*((*event_it)->getObjectIt()))->getCreationTime());
    }
    // Delete polaron and its events
    deleteObject((*event_it)->getObjectIt());
    // Update polaron counters
    if(!charge){
        N_electrons_collected++;
        N_electrons--;
    }
    else{
        N_holes_collected++;
        N_holes--;
    }
    // Log event
    if(loggingEnabled()){
        ostringstream msg;
        if(!charge){
            msg << "Electron " << polaron_tag << " was extracted from site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
        }
        else{
            msg << "Hole " << polaron_tag << " was extracted from site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
        }
        logMSG(msg);
    }
    // Find all nearby objects and calculate their events
    auto neighbors = findRecalcNeighbors(coords_initial);
    calculateObjectListEvents(neighbors);
    return true;
}

bool OSC_Sim::executePolaronHop(const list<unique_ptr<Event>>::iterator event_it){
    if(isOccupied((*event_it)->getDestCoords())){
        cout << "Polaron hop cannot be executed. Destination site is already occupied." << endl;
        return false;
    }
    else{
        // Log event
        if(loggingEnabled()){
            ostringstream msg;
            msg << "Polaron " << (*((*event_it)->getObjectIt()))->getTag() << " hopping to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
            logMSG(msg);
        }
        return executeObjectHop(event_it);
    }
}

bool OSC_Sim::executePolaronRecombination(const list<unique_ptr<Event>>::iterator event_it){
    // Get event info
    int polaron_tag = (*((*event_it)->getObjectIt()))->getTag();
    int target_tag = (*((*event_it)->getObjectTargetIt()))->getTag();
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    Coords coords_dest = (*event_it)->getDestCoords();
    // Update simulation time
    incrementTime((*event_it)->getWaitTime());
    // Delete polaron and its events
    deleteObject((*event_it)->getObjectIt());
    deleteObject((*event_it)->getObjectTargetIt());
    // Update polaron counters
    N_electrons_recombined++;
    N_holes_recombined++;
    N_electrons--;
    N_holes--;
    if(polaron_tag==target_tag){
        N_geminate_recombinations++;
    }
    else{
        N_bimolecular_recombinations++;
    }
    // Log event
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Electron " << polaron_tag << " at site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << " recombined with hole " << target_tag << " at site " << coords_dest.x << "," << coords_dest.y << "," << coords_dest.z << "." << endl;
        logMSG(msg);
    }
    // Find all nearby objects
    auto neighbors = findRecalcNeighbors(coords_initial);
    auto neighbors2 = findRecalcNeighbors(coords_dest);
    neighbors.insert(neighbors.end(),neighbors2.begin(),neighbors2.end());
    removeObjectItDuplicates(neighbors);
    // Calculate events for all nearby objects
    calculateObjectListEvents(neighbors);
    return true;
}

void OSC_Sim::generateElectron(const Coords& coords,int tag=0){
    if(tag==0){
        tag = N_electrons_created+1;
    }
    // Create the new electron and add it to the simulation
    Polaron electron_new(getTime(),tag,coords,false);
    electrons.push_back(electron_new);
    auto object_it = addObject(&electrons.back());
    // Add placeholder events to the corresponding lists
    Polaron_Hop hop_event;
    electron_hop_events.push_back(hop_event);
    Polaron_Recombination recombination_event;
    recombination_event.setObjectIt(object_it);
    polaron_recombination_events.push_back(recombination_event);
    Polaron_Extraction extraction_event;
    extraction_event.setObjectIt(object_it);
    electron_extraction_events.push_back(extraction_event);
    // Update exciton counters
    N_electrons_created++;
    N_electrons++;
    // Log event
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Created electron " << electron_new.getTag() << " at site " << coords.x << "," << coords.y << "," << coords.z << "." << endl;
        logMSG(msg);
    }
}

void OSC_Sim::generateHole(const Coords& coords,int tag=0){
    if(tag==0){
        tag = N_holes_created+1;
    }
    // Create the new hole and add it to the simulation
    Polaron hole_new(getTime(),tag,coords,true);
    holes.push_back(hole_new);
    auto object_it = addObject(&holes.back());
    // Add placeholder events to the corresponding lists
    Polaron_Hop hop_event;
    hole_hop_events.push_back(hop_event);
    Polaron_Extraction extraction_event;
    extraction_event.setObjectIt(object_it);
    hole_extraction_events.push_back(extraction_event);
    // Update exciton counters
    N_electrons_created++;
    N_electrons++;
    // Log event
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Created hole " << hole_new.getTag() << " at site " << coords.x << "," << coords.y << "," << coords.z << "." << endl;
        logMSG(msg);
    }
}

void OSC_Sim::generateToFPolarons(){
    int num = 0;
    Coords coords;
    // Create electrons at the top plane of the lattice
    if(!Polaron_type){
        coords.z = getHeight()-1;
    }
    // Create holes at the bottom plane of the lattice
    else{
        coords.z = 0;
    }
    while(num<ToF_initial_polarons){
        coords.x = getRandomX();
        coords.y = getRandomY();
        // If the site is already occupied, pick a new site
        if(isOccupied(coords)){
            continue;
        }
        // If phase restriction is enabled, electrons cannot be created on donor sites
        else if(Enable_phase_restriction && !Polaron_type && getSiteType(coords)==(short)1){
            continue;
        }
        // If phase restriction is enabled, holes cannot be created on acceptor sites
        else if(Enable_phase_restriction && Polaron_type && getSiteType(coords)==(short)2){
            continue;
        }
        else if(!Polaron_type){
            generateElectron(coords);
            num++;
        }
        else{
            generateHole(coords);
            num++;
        }
    }
    auto object_its = getAllObjectIts();
    calculateObjectListEvents(object_its);
}

vector<double> OSC_Sim::getDiffusionData(){
    return diffusion_distances;
}

list<Exciton>::iterator OSC_Sim::getExcitonIt(const unique_ptr<Object>& object_ptr){
    for(auto it=excitons.begin();it!=excitons.end();++it){
        if(object_ptr->getTag()==it->getTag()){
            return it;
        }
    }
}

int OSC_Sim::getN_excitons_created(){
    return N_excitons_created;
}

int OSC_Sim::getN_excitons_recombined(){
    return N_excitons_recombined;
}

bool OSC_Sim::getObjectCharge(const list<unique_ptr<Object>>::iterator object_it){
    if((*object_it)->getName().compare(Polaron::name)==0){
        return static_cast<Polaron*>(object_it->get())->getCharge();
    }
}

list<Polaron>::iterator OSC_Sim::getPolaronIt(const unique_ptr<Object>& object_ptr){
    if(object_ptr->getName().compare(Polaron::name)==0){
        // electrons
        if(!static_cast<Polaron*>(object_ptr.get())->getCharge()){
            for(auto it=electrons.begin();it!=electrons.end();++it){
                if(object_ptr->getTag()==it->getTag()){
                    return it;
                }
            }
        }
        // holes
        else{
            for(auto it=holes.begin();it!=holes.end();++it){
                if(object_ptr->getTag()==it->getTag()){
                    return it;
                }
            }
        }
    }
}

double OSC_Sim::getSiteEnergy(const Coords& coords){
    return sites[getSiteIndex(coords)].getEnergy();
}

short OSC_Sim::getSiteType(const Coords& coords){
    return static_cast<Site_OSC*>(getSiteIt(coords)->get())->getType();
}

void OSC_Sim::initializeArchitecture(){
    N_donor_sites = 0;
    N_acceptor_sites = 0;
    if(Enable_neat){
        N_donor_sites = getNumSites();
        N_acceptor_sites = 0;
        for(auto it=sites.begin();it!=sites.end();++it){
            it->setType(1);
        }
    }
    else if(Enable_bilayer){
        Coords coords;
        for(int x=0;x<getLength();x++){
            for(int y=0;y<getWidth();y++){
                for(int z=0;z<getHeight();z++){
                    coords.x = x;
                    coords.y = y;
                    coords.z = z;
                    if(z<Thickness_acceptor){
                        sites[getSiteIndex(coords)].setType(2);
                        N_acceptor_sites++;
                    }
                    else{
                        sites[getSiteIndex(coords)].setType(1);
                        N_donor_sites++;
                    }
                }
            }
        }
    }
    else if(Enable_random_blend){
        vector<short> site_types;
        site_types.assign(getNumSites(),1);
        for(int i=0;i<(int)getNumSites()*Acceptor_conc;i++){
            site_types[i] = 2;
            N_acceptor_sites++;
        }
        N_donor_sites = getNumSites()-N_acceptor_sites;
        shuffle(site_types.begin(),site_types.end(),gen);
    }
    if(Enable_gaussian_dos){
        site_energies_donor.assign(N_donor_sites,0);
        site_energies_acceptor.assign(N_acceptor_sites,0);
        createGaussianDOSVector(site_energies_donor,0,Energy_stdev_donor);
        createGaussianDOSVector(site_energies_acceptor,0,Energy_stdev_acceptor);
    }
    else if(Enable_exponential_dos){
        site_energies_donor.assign(N_donor_sites,0);
        site_energies_donor.assign(N_acceptor_sites,0);
        createExponentialDOSVector(site_energies_donor,0,Energy_urbach_donor);
        createExponentialDOSVector(site_energies_acceptor,0,Energy_urbach_acceptor);
    }
    else{
        site_energies_donor.push_back(0);
        site_energies_acceptor.push_back(0);
    }
    unique_ptr<Site> site_ptr;
    int donor_count = 0;
    int acceptor_count = 0;
    for(int i=0;i<getNumSites();i++){
        if(Enable_gaussian_dos || Enable_exponential_dos){
            if(sites[i].getType()==(short)1){
                sites[i].setEnergyIt(site_energies_donor.begin()+donor_count);
                donor_count++;
            }
            else if(sites[i].getType()==(short)2){
                sites[i].setEnergyIt(site_energies_acceptor.begin()+acceptor_count);
                acceptor_count++;
            }
            else{
                cout << "Error! Undefined site type detected while assigning site energies." << endl;
            }
        }
        else{
            if(sites[i].getType()==(short)1){
                sites[i].setEnergyIt(site_energies_donor.begin());
            }
            else if(sites[i].getType()==(short)2){
                sites[i].setEnergyIt(site_energies_acceptor.begin());
            }
            else{
                cout << "Error! Undefined site type detected while assigning site energies." << endl;
            }
        }
        addSite(&sites[i]);
    }
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

bool OSC_Sim::siteContainsHole(const Coords& coords){
    auto object_it = (*getSiteIt(coords))->getObjectIt();
    if((*object_it)->getName().compare(Polaron::name)==0){
        return static_cast<Polaron*>(object_it->get())->getCharge();
    }
    return false;
}
