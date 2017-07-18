// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "OSC_Sim.h"

OSC_Sim::OSC_Sim() {

}

OSC_Sim::~OSC_Sim() {

}

void OSC_Sim::init(const Parameters_OPV& params,const int id){
    // Set parameters of Simulation base class
    Simulation::init(params,id);
    // Set Additional General Parameters
    Bias = params.Bias;
    // Morphology Parameters
    Enable_neat = params.Enable_neat;
    Enable_bilayer = params.Enable_bilayer;
    Thickness_donor = params.Thickness_donor;
    Thickness_acceptor = params.Thickness_acceptor;
    Enable_random_blend = params.Enable_random_blend;
    Acceptor_conc = params.Acceptor_conc;
    Enable_import_morphology = params.Enable_import_morphology;
    Morphology_file = params.Morphology_file;
    // Test Parameters
    N_tests = params.N_tests;
    Enable_exciton_diffusion_test = params.Enable_exciton_diffusion_test;
    Enable_ToF_test = params.Enable_ToF_test;
    ToF_polaron_type = params.ToF_polaron_type;
    ToF_initial_polarons = params.ToF_initial_polarons;
    ToF_transient_start = params.ToF_transient_start;
    ToF_transient_end = params.ToF_transient_end;
    ToF_pnts_per_decade = params.ToF_pnts_per_decade;
    Enable_IQE_test = params.Enable_IQE_test;
    IQE_time_cutoff = params.IQE_time_cutoff;
    Enable_dynamics_test = params.Enable_dynamics_test;
    Enable_dynamics_extraction = params.Enable_dynamics_extraction;
    Dynamics_initial_exciton_conc = params.Dynamics_initial_exciton_conc;
    Dynamics_transient_start = params.Dynamics_transient_start;
    Dynamics_transient_end = params.Dynamics_transient_end;
    Dynamics_pnts_per_decade = params.Dynamics_pnts_per_decade;
    // Exciton Parameters
    Exciton_generation_rate_donor = params.Exciton_generation_rate_donor;
    Exciton_generation_rate_acceptor = params.Exciton_generation_rate_acceptor;
    Exciton_lifetime_donor = params.Exciton_lifetime_donor;
    Exciton_lifetime_acceptor = params.Exciton_lifetime_acceptor;
    R_exciton_hopping_donor = params.R_exciton_hopping_donor;
    R_exciton_hopping_acceptor = params.R_exciton_hopping_acceptor;
    FRET_cutoff = params.FRET_cutoff;
    E_exciton_binding_donor = params.E_exciton_binding_donor;
    E_exciton_binding_acceptor = params.E_exciton_binding_acceptor;
    R_exciton_dissociation_donor = params.R_exciton_dissociation_donor;
    R_exciton_dissociation_acceptor = params.R_exciton_dissociation_acceptor;
    Exciton_dissociation_cutoff = params.Exciton_dissociation_cutoff;
    // Polaron Parameters
    Enable_phase_restriction = params.Enable_phase_restriction;
    R_polaron_hopping_donor = params.R_polaron_hopping_donor;
    R_polaron_hopping_acceptor = params.R_polaron_hopping_acceptor;
    Polaron_localization_donor = params.Polaron_localization_donor;
    Polaron_localization_acceptor = params.Polaron_localization_acceptor;
    Enable_miller_abrahams = params.Enable_miller_abrahams;
    Enable_marcus = params.Enable_marcus;
    Reorganization_donor = params.Reorganization_donor;
    Reorganization_acceptor = params.Reorganization_acceptor;
    R_polaron_recombination = params.R_polaron_recombination;
    Polaron_hopping_cutoff = params.Polaron_hopping_cutoff;
    Enable_gaussian_polaron_delocalization = params.Enable_gaussian_polaron_delocalization;
    Polaron_delocalization_length = params.Polaron_delocalization_length;
    // Additional Lattice Parameters
    Homo_donor = params.Homo_donor;
    Lumo_donor = params.Lumo_donor;
    Homo_acceptor = params.Homo_acceptor;
    Lumo_acceptor = params.Lumo_acceptor;
    Enable_gaussian_dos = params.Enable_gaussian_dos;
    Energy_stdev_donor = params.Energy_stdev_donor;
    Energy_stdev_acceptor = params.Energy_stdev_acceptor;
    Enable_exponential_dos = params.Enable_exponential_dos;
    Energy_urbach_donor = params.Energy_urbach_donor;
    Energy_urbach_acceptor = params.Energy_urbach_acceptor;
    // Coulomb Calculation Parameters
    Dielectric_donor = params.Dielectric_donor;
    Dielectric_acceptor = params.Dielectric_acceptor;
    Coulomb_cutoff = params.Coulomb_cutoff;
    // Output files
    //
    // Additional Parameters
    Error_found = false;
    // Initialize counters
    N_excitons = 0;
    N_excitons_created = 0;
    N_excitons_created_donor = 0;
    N_excitons_created_acceptor = 0;
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
    site.setType(0);
    sites.assign(lattice.getNumSites(),site);
    // Initialize Film Architecture and Site Energies
    initializeArchitecture();
    // Calculate Coulomb interactions lookup table
    double avgDielectric = (Dielectric_donor+Dielectric_acceptor)/2;
    int range = ceil((Coulomb_cutoff/lattice.getUnitSize())*(Coulomb_cutoff/lattice.getUnitSize()));
    Coulomb_table.assign(range+1,0);
    for(int i=1;i<(int)Coulomb_table.size();i++){
        Coulomb_table[i] = ((Coulomb_constant*Elementary_charge)/avgDielectric)/(1e-9*lattice.getUnitSize()*sqrt((double)i));
        if(Enable_gaussian_polaron_delocalization){
            Coulomb_table[i] *= erf((lattice.getUnitSize()*sqrt((double)i))/(Polaron_delocalization_length*sqrt(2)));
        }
    }
    // Initialize electrical potential vector
    E_potential.assign(lattice.getHeight(),0);
    for(int i=0;i<lattice.getHeight();i++){
        E_potential[i] = (Bias*lattice.getHeight()/(lattice.getHeight()+1))-(Bias/(lattice.getHeight()+1))*i;
    }
    // Initialize exciton creation event
    R_exciton_generation_donor = Exciton_generation_rate_donor*N_donor_sites*intpow(1e-7*lattice.getUnitSize(),3);
    R_exciton_generation_acceptor = Exciton_generation_rate_acceptor*N_acceptor_sites*intpow(1e-7*lattice.getUnitSize(),3);
    if(Enable_exciton_diffusion_test || Enable_IQE_test){
        isLightOn = true;
        exciton_creation_event.calculateExecutionTime(R_exciton_generation_donor+R_exciton_generation_acceptor,getTime());
        exciton_creation_it = addEvent(&exciton_creation_event);
    }
    else if(Enable_dynamics_test){
        // Initialize data structures
        isLightOn = false;
        double step_size = 1.0/(double)Dynamics_pnts_per_decade;
        int num_steps = floor((log10(Dynamics_transient_end)-log10(Dynamics_transient_start))/step_size);
        transient_times.assign(num_steps,0);
        for(int i=0;i<(int)transient_times.size();i++){
            transient_times[i] = pow(10,log10(Dynamics_transient_start)+i*step_size);
        }
        transient_excitons.assign(num_steps,0);
        transient_electrons.assign(num_steps,0);
        transient_holes.assign(num_steps,0);
        // Create initial test excitons
        generateDynamicsExcitons();
    }
    else if(Enable_ToF_test){
        // Initialize data structures
        isLightOn = false;
        double step_size = 1.0/(double)ToF_pnts_per_decade;
        int num_steps = floor((log10(ToF_transient_end)-log10(ToF_transient_start))/step_size);
        transient_times.assign(num_steps,0);
        for(int i=0;i<(int)transient_times.size();i++){
            transient_times[i] = pow(10,log10(ToF_transient_start)+(i+0.5)*step_size);
        }
        transient_counts.assign(num_steps,0);
        transient_velocities.assign(num_steps,0);
        transient_energies.assign(num_steps,0);
        // Create initial test polarons
        generateToFPolarons();
    }
}

double OSC_Sim::calculateCoulomb(const list<Object*>::iterator object_it,const Coords& coords){
    bool charge;
    if((*object_it)->getName().compare(Polaron::name)==0){
        charge = static_cast<Polaron*>(*object_it)->getCharge();
    }
    else{
        cout << getId() << ": Error! Coulomb interactions calculated for an object that is not a Polaron." << endl;
        Error_found = true;
    }
    static const double avgDielectric = (Dielectric_donor+Dielectric_acceptor)/2;
    static const double image_interactions = (Elementary_charge/(16*Pi*avgDielectric*Vacuum_permittivity))*1e9;
    double Energy = 0;
    double distance;
    int distance_sq_lat;
    static const int range = ceil((Coulomb_cutoff/ lattice.getUnitSize())*(Coulomb_cutoff/ lattice.getUnitSize()));
    // Loop through electrons
    for(auto it=electrons.begin();it!=electrons.end();++it){
        distance_sq_lat = lattice.calculateLatticeDistanceSquared(coords,it->getCoords());
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
        distance_sq_lat = lattice.calculateLatticeDistanceSquared(coords,it->getCoords());
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
    if(!lattice.isZPeriodic()){
        distance = lattice.getUnitSize()*((double)(lattice.getHeight()-coords.z)-0.5);
        if(!((distance-0.0001)>Coulomb_cutoff)){
            Energy -= image_interactions/distance;
        }
        distance = lattice.getUnitSize()*((double)(coords.z+1)-0.5);
        if(!((distance-0.0001)>Coulomb_cutoff)){
            Energy -= image_interactions/distance;
        }
    }
    return Energy;
}

double OSC_Sim::calculateCoulomb(const bool charge,const Coords& coords){
    static const double avgDielectric = (Dielectric_donor+Dielectric_acceptor)/2;
    static const double image_interactions = (Elementary_charge/(16*Pi*avgDielectric*Vacuum_permittivity))*1e9;
    double Energy = 0;
    double distance;
    int distance_sq_lat;
    static const int range = ceil((Coulomb_cutoff/ lattice.getUnitSize())*(Coulomb_cutoff/ lattice.getUnitSize()));
    // Loop through electrons
    for(auto it=electrons.begin();it!=electrons.end();++it){
        distance_sq_lat = lattice.calculateLatticeDistanceSquared(coords,it->getCoords());
        if(!distance_sq_lat>range){
            if(!charge){
                Energy += Coulomb_table[distance_sq_lat];
            }
            else{
                Energy -= Coulomb_table[distance_sq_lat];
            }
        }
    }
    // Loop through holes
    for(auto it=holes.begin();it!=holes.end();++it){
        distance_sq_lat = lattice.calculateLatticeDistanceSquared(coords,it->getCoords());
        if(!distance_sq_lat>range){
            if(charge){
                Energy += Coulomb_table[distance_sq_lat];
            }
            else{
                Energy -= Coulomb_table[distance_sq_lat];
            }
        }
    }
    // Add electrode image charge interactions
    if(!lattice.isZPeriodic()){
        distance = lattice.getUnitSize()*((double)(lattice.getHeight()-coords.z)-0.5);
        if(!((distance-0.0001)>Coulomb_cutoff)){
            Energy -= image_interactions/distance;
        }
        distance = lattice.getUnitSize()*((double)(coords.z+1)-0.5);
        if(!((distance-0.0001)>Coulomb_cutoff)){
            Energy -= image_interactions/distance;
        }
    }
    return Energy;
}

double OSC_Sim::calculateDiffusionLength_avg(){
    return vector_avg(diffusion_distances);
}

double OSC_Sim::calculateDiffusionLength_stdev(){
    return vector_stdev(diffusion_distances);
}

double OSC_Sim::calculateMobility_avg(){
    auto mobilities = transit_times;
    for(int i=0;i<(int)mobilities.size();i++){
        mobilities[i] = intpow(1e-7*lattice.getUnitSize()*lattice.getHeight(),2)/(fabs(Bias)*transit_times[i]);
    }
    return vector_avg(mobilities);
}

double OSC_Sim::calculateMobility_stdev(){
    auto mobilities = transit_times;
    for(int i=0;i<(int)mobilities.size();i++){
        mobilities[i] = intpow(1e-7*lattice.getUnitSize()*lattice.getHeight(),2)/(fabs(Bias)*transit_times[i]);
    }
    return vector_stdev(mobilities);
}

vector<double> OSC_Sim::calculateTransitTimeDist(const vector<double>& data,const int counts){
    double step_size = 1.0/(double)ToF_pnts_per_decade;
    vector<double> dist(transient_times.size(),0);
    for(int i=0;i<(int)data.size();i++){
        for(int j=0;j<(int)transient_times.size();j++){
            if(data[i] > pow(10,log10(transient_times[j])-0.5*step_size) && data[i] < pow(10,log10(transient_times[j])+0.5*step_size)){
                dist[j] += 1.0/counts;
            }
        }
    }
    return dist;
}

double OSC_Sim::calculateTransitTime_avg(){
    return vector_avg(transit_times);
}

double OSC_Sim::calculateTransitTime_stdev(){
    return vector_stdev(transit_times);
}

Coords OSC_Sim::calculateExcitonCreationCoords(){
    static uniform_real_distribution<double> dist(0.0,R_exciton_generation_donor+R_exciton_generation_acceptor);
    double num = dist(gen);
    short type_target;
    if(num<R_exciton_generation_donor){
        type_target = 1;
    }
    else{
        type_target = 2;
    }
    Coords dest_coords;
    int N_tries = 0;
    while(N_tries<1000000){
        dest_coords = lattice.getRandomCoords();
        if(loggingEnabled()){
            *Logfile << "Attempting to create exciton at " << dest_coords.x << "," << dest_coords.y << "," << dest_coords.z << "." << endl;
        }
        N_tries++;
        if(!lattice.isOccupied(dest_coords) && getSiteType(dest_coords)==type_target){
            return dest_coords;
        }
    }
    cout << getId() << ": Error! An empty site for exciton creation could not be found." << endl;
    Error_found = true;
    return dest_coords;
}

void OSC_Sim::calculateExcitonEvents(const list<Object*>::iterator object_it){
    const auto exciton_it = getExcitonIt(*object_it);
    const Coords object_coords = exciton_it->getCoords();
    if(loggingEnabled()){
        *Logfile << "Calculating events for exciton " << exciton_it->getTag() << " at site " << object_coords.x << "," << object_coords.y << "," << object_coords.z << "." << endl;
    }
    Coords dest_coords;
    double distance,E_delta,Coulomb_final,rate;
    int index;
    // Exciton hopping and dissociation
    static const int range = ceil( ((FRET_cutoff>Exciton_dissociation_cutoff) ? (FRET_cutoff):(Exciton_dissociation_cutoff))/ lattice.getUnitSize());
    static const int dim = (2*range+1);
    static vector<Exciton_Hop> hops_temp(dim*dim*dim);
    static vector<Exciton_Dissociation> dissociations_temp(dim*dim*dim);
    static vector<bool> hops_valid(dim*dim*dim,false);
    static vector<bool> dissociations_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                if(!checkMoveEventValidity(object_coords,i,j,k)){
                    hops_valid[index] = false;
                    dissociations_valid[index] = false;
                    continue;
                }
                dest_coords = lattice.calculateDestinationCoords(object_coords,i,j,k);
                if(lattice.isOccupied(dest_coords)){
                    hops_valid[index] = false;
                    dissociations_valid[index] = false;
                    continue;
                }
                distance = lattice.getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                // Dissociation event
                if(getSiteType(object_coords)!=getSiteType(dest_coords) && !((distance-0.0001)>Exciton_dissociation_cutoff)){
                    dissociations_temp[index].setObjectIt(object_it);
                    dissociations_temp[index].setDestCoords(dest_coords);
                    // Exciton is starting from a donor site
                    if(getSiteType(object_coords)==(short)1){
                        Coulomb_final = calculateCoulomb(true,object_coords)+calculateCoulomb(false,dest_coords)-Coulomb_table[i*i+j*j+k*k];
                        E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords))-(Lumo_acceptor-Lumo_donor)+(Coulomb_final+E_exciton_binding_donor)+(E_potential[dest_coords.z]-E_potential[object_coords.z]);
                        if(Enable_miller_abrahams){
                            dissociations_temp[index].calculateExecutionTime(R_exciton_dissociation_donor,Polaron_localization_donor,distance,E_delta,getTemperature(),getTime());
                        }
                        else{
                            dissociations_temp[index].calculateExecutionTime(R_exciton_dissociation_donor,Polaron_localization_donor,distance,E_delta,Reorganization_donor,getTemperature(),getTime());
                        }
                    }
                    // Exciton is starting from an acceptor site
                    else{
                        Coulomb_final = calculateCoulomb(false,object_coords)+calculateCoulomb(true,dest_coords)-Coulomb_table[i*i+j*j+k*k];
                        E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords))+(Homo_donor-Homo_acceptor)+(Coulomb_final+E_exciton_binding_donor)-(E_potential[dest_coords.z]-E_potential[object_coords.z]);
                        if(Enable_miller_abrahams){
                            dissociations_temp[index].calculateExecutionTime(R_exciton_dissociation_acceptor,Polaron_localization_acceptor,distance,E_delta,getTemperature(),getTime());
                        }
                        else{
                            dissociations_temp[index].calculateExecutionTime(R_exciton_dissociation_acceptor,Polaron_localization_acceptor,distance,E_delta,Reorganization_acceptor,getTemperature(),getTime());
                        }
                    }
                    dissociations_valid[index] = true;
                }
                else{
                    dissociations_valid[index] = false;
                }
                // Hop event
                if(getSiteType(object_coords)==getSiteType(dest_coords) && !((distance-0.0001)>FRET_cutoff)){
                    hops_temp[index].setObjectIt(object_it);
                    hops_temp[index].setDestCoords(dest_coords);
                    E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords));
                    if(getSiteType(object_coords)==(short)1){
                        hops_temp[index].calculateExecutionTime(R_exciton_hopping_donor,distance,E_delta,getTemperature(),getTime());
                    }
                    else{
                        hops_temp[index].calculateExecutionTime(R_exciton_hopping_acceptor,distance,E_delta,getTemperature(),getTime());
                    }
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
    if(getSiteType(object_coords)==(short)1){
        rate = 1/Exciton_lifetime_donor;
    }
    else if(getSiteType(object_coords)==(short)2){
        rate = 1/Exciton_lifetime_acceptor;
    }
    recombination_event_it->calculateExecutionTime(rate,getTime());
    // Determine the fastest valid hop event
    bool No_hops_valid = true;
    auto hop_target_it = hops_temp.end();
    for(auto it=hops_temp.begin();it!=hops_temp.end();++it){
        if(hops_valid[std::distance(hops_temp.begin(),it)] && (hop_target_it==hops_temp.end() || it->getExecutionTime()<hop_target_it->getExecutionTime())){
            hop_target_it = it;
            No_hops_valid = false;
        }
    }
    // Determine the fastest valid dissociation event
    bool No_dissociations_valid = true;
    auto dissociation_target_it = dissociations_temp.end();
    for(auto it=dissociations_temp.begin();it!=dissociations_temp.end();++it){
        if(dissociations_valid[std::distance(dissociations_temp.begin(),it)] && (dissociation_target_it==dissociations_temp.end() || it->getExecutionTime()<dissociation_target_it->getExecutionTime())){
            dissociation_target_it = it;
            No_dissociations_valid = false;
        }
    }
    // Determine fastest valid event
    int selection = 1;
    double best_time = recombination_event_it->getExecutionTime();
    if(!No_hops_valid && hop_target_it->getExecutionTime()<best_time){
        selection = 2;
        best_time = hop_target_it->getExecutionTime();
    }
    if(!No_dissociations_valid && dissociation_target_it->getExecutionTime()<best_time){
        selection = 3;
    }
    switch(selection){
        // Recombination event is fastest
        case 1:{
            setEvent(exciton_it->getEventIt(),&(*recombination_event_it));
            break;
        }
        // Hop event is fastest
        case 2:{
            auto hop_list_it = exciton_hop_events.begin();
            advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
            *hop_list_it = *hop_target_it;
            setEvent(exciton_it->getEventIt(),&(*hop_list_it));
            break;
        }
        // Dissociation event is fastest
        case 3:{
            auto dissociation_list_it = exciton_dissociation_events.begin();
            advance(dissociation_list_it,std::distance(excitons.begin(),exciton_it));
            *dissociation_list_it = *dissociation_target_it;
            setEvent(exciton_it->getEventIt(),&(*dissociation_list_it));
            break;
        }
        // No valid event found
        default:{
            cout << getId() << ": Error! No valid events could be calculated." << endl;
            Error_found = true;
            return;
        }
    }
}

void OSC_Sim::calculateObjectListEvents(const vector<list<Object*>::iterator>& object_it_vec){
    if(loggingEnabled()){
        *Logfile << "Calculating events for " << object_it_vec.size() << " objects:" << endl;
    }
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

void OSC_Sim::calculatePolaronEvents(const list<Object*>::iterator object_it){
    const auto polaron_it = getPolaronIt(*object_it);
    const Coords object_coords = polaron_it->getCoords();
    if(loggingEnabled()){
        if(!polaron_it->getCharge()){
            *Logfile << "Calculating events for electron " << polaron_it->getTag() << " at site " << object_coords.x << "," << object_coords.y << "," << object_coords.z << "." << endl;
        }
        else{
            *Logfile << "Calculating events for hole " << polaron_it->getTag() << " at site " << object_coords.x << "," << object_coords.y << "," << object_coords.z << "." << endl;
        }
    }
    Coords dest_coords;
    double distance,E_delta;
    int index;
    double Coulomb_i = calculateCoulomb(object_it,object_coords);
    // Calculate Polaron hopping and recombination events
    static const int range = ceil(Polaron_hopping_cutoff/ lattice.getUnitSize());
    static const int dim = (2*range+1);
    static vector<Polaron_Hop> hops_temp(dim*dim*dim);
    static vector<Polaron_Recombination> recombinations_temp(dim*dim*dim);
    static vector<bool> hops_valid(dim*dim*dim,false);
    static vector<bool> recombinations_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                hops_valid[index] = false;
                recombinations_valid[index] = false;
                if(!checkMoveEventValidity(object_coords,i,j,k)){
                    continue;
                }
                dest_coords = lattice.calculateDestinationCoords(object_coords,i,j,k);
                // Recombination events
                // If destination site is occupied by a hole Polaron and the main Polaron is an electron, check for a possible recombination event
                if(lattice.isOccupied(dest_coords) && !polaron_it->getCharge() && siteContainsHole(dest_coords)){
                    distance = lattice.getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                    if(!((distance-0.0001)>Polaron_hopping_cutoff)){
                        recombinations_temp[index].setObjectIt(object_it);
                        if(getSiteType(object_coords)==(short)1){
                            recombinations_temp[index].calculateExecutionTime(R_polaron_recombination,Polaron_localization_donor,distance,0,getTemperature(),getTime());
                        }
                        else if(getSiteType(object_coords)==(short)2){
                            recombinations_temp[index].calculateExecutionTime(R_polaron_recombination,Polaron_localization_acceptor,distance,0,getTemperature(),getTime());
                        }
                        recombinations_temp[index].setDestCoords(dest_coords);
                        recombinations_temp[index].setObjectTargetIt((*lattice.getSiteIt(dest_coords))->getObjectIt());
                        recombinations_valid[index] = true;
                    }
                }
                // Hop events
                // If destination site is unoccupied and either phase restriction is disabled or the starting site and destination sites have the same type, check for a possible hop event
                else if(!lattice.isOccupied(dest_coords) && (!Enable_phase_restriction || getSiteType(object_coords)==getSiteType(dest_coords))){
                    distance = lattice.getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                    if(!((distance-0.0001)>Polaron_hopping_cutoff)){
                        hops_temp[index].setObjectIt(object_it);
                        hops_temp[index].setDestCoords(dest_coords);
                        E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords))+(calculateCoulomb(object_it,dest_coords)-Coulomb_i);
                        if(!getObjectCharge(object_it)){
                            E_delta += (E_potential[dest_coords.z]-E_potential[object_coords.z]);
                        }
                        else{
                            E_delta -= (E_potential[dest_coords.z]-E_potential[object_coords.z]);
                        }
                        if(getSiteType(object_coords)==(short)1){
                            if(getSiteType(dest_coords)==(short)2){
                                if(!getObjectCharge(object_it)){
                                    E_delta -= (Lumo_acceptor-Lumo_donor);
                                }
                                else{
                                    E_delta -= (Homo_acceptor-Homo_donor);
                                }
                            }
                            if(Enable_miller_abrahams){
                                hops_temp[index].calculateExecutionTime(R_polaron_hopping_donor,Polaron_localization_donor,distance,E_delta,getTemperature(),getTime());
                            }
                            else{
                                hops_temp[index].calculateExecutionTime(R_polaron_hopping_donor,Polaron_localization_donor,distance,E_delta,Reorganization_donor,getTemperature(),getTime());
                            }
                        }
                        else if(getSiteType(object_coords)==(short)2){
                            if(getSiteType(dest_coords)==(short)1){
                                if(!getObjectCharge(object_it)){
                                    E_delta -= (Lumo_donor-Lumo_acceptor);
                                }
                                else{
                                    E_delta -= (Homo_donor-Homo_acceptor);
                                }
                            }
                            if(Enable_miller_abrahams){
                                hops_temp[index].calculateExecutionTime(R_polaron_hopping_acceptor,Polaron_localization_acceptor,distance,E_delta,getTemperature(),getTime());
                            }
                            else{
                                hops_temp[index].calculateExecutionTime(R_polaron_hopping_acceptor,Polaron_localization_acceptor,distance,E_delta,Reorganization_acceptor,getTemperature(),getTime());
                            }
                        }
                        hops_valid[index] = true;
                    }
                }
            }
        }
    }
    // Calculate possible extraction event
    // Electrons are extracted at the bottom of the lattice (z=-1)
    // Holes are extracted at the top of the lattice (z=Height)
    bool No_extraction_valid = true;
    list<Polaron_Extraction>::iterator extraction_event_it;
    if(!Enable_dynamics_test || Enable_dynamics_extraction){
        // If electron, charge is false
        if(!polaron_it->getCharge()){
            distance = lattice.getUnitSize()*((double)(object_coords.z+1)-0.5);
            if(!((distance-0.0001)>Polaron_hopping_cutoff)){
                extraction_event_it = electron_extraction_events.begin();
                advance(extraction_event_it,std::distance(electrons.begin(),polaron_it));
                No_extraction_valid = false;
            }
        }
        // If hole, charge is true
        else{
            distance = lattice.getUnitSize()*((double)(lattice.getHeight()-object_coords.z)-0.5);
            if(!((distance-0.0001)>Polaron_hopping_cutoff)){
                extraction_event_it = hole_extraction_events.begin();
                advance(extraction_event_it,std::distance(holes.begin(),polaron_it));
                No_extraction_valid = false;
            }
        }
        if(!No_extraction_valid){
            if(getSiteType(object_coords)==(short)1){
                extraction_event_it->calculateExecutionTime(R_polaron_hopping_donor,distance,Polaron_localization_donor,0,getTemperature(),getTime());
            }
            else if(getSiteType(object_coords)==(short)2){
                extraction_event_it->calculateExecutionTime(R_polaron_hopping_acceptor,distance,Polaron_localization_acceptor,0,getTemperature(),getTime());
            }
        }
    }
    // Determine the fastest hop event
    bool No_hops_valid = true;
    auto hop_target_it = hops_temp.end();
    for(auto it=hops_temp.begin();it!=hops_temp.end();++it){
        if(hops_valid[std::distance(hops_temp.begin(),it)] && (hop_target_it==hops_temp.end() || it->getExecutionTime()<hop_target_it->getExecutionTime())){
            hop_target_it = it;
            No_hops_valid = false;
        }
    }
    // Determine the fastest recombination event
    bool No_recombinations_valid = true;
    auto recombination_target_it = recombinations_temp.end();
    for(auto it=recombinations_temp.begin();it!=recombinations_temp.end();++it){
        if(recombinations_valid[std::distance(recombinations_temp.begin(),it)] && (recombination_target_it==recombinations_temp.end() || it->getExecutionTime()<recombination_target_it->getExecutionTime())){
            recombination_target_it = it;
            No_recombinations_valid = false;
        }
    }
    // If no valid event is found, return
    if(No_extraction_valid && No_hops_valid && No_recombinations_valid){
        setEvent(polaron_it->getEventIt(),nullptr);
        return;
    }
    // Determine the fastest valid event
    short selection = -1;
    double best_time = -1;
    if(!No_extraction_valid){
        selection = 1;
        best_time = extraction_event_it->getExecutionTime();
    }
    if(!No_hops_valid && (best_time<0 || hop_target_it->getExecutionTime()<best_time)){
        selection = 2;
        best_time = hop_target_it->getExecutionTime();
    }
    if(!No_recombinations_valid && (best_time<0 || recombination_target_it->getExecutionTime()<best_time)){
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
                cout << getId() << ": Error! Only electrons can initiate polaron recombination." << endl;
                Error_found = true;
                return;
            }
            *recombination_list_it = *recombination_target_it;
            setEvent(polaron_it->getEventIt(),&(*recombination_list_it));
            break;
        }
        // No valid event found
        default:{
            cout << getId() << ": Error! No valid events could be calculated." << endl;
            Error_found = true;
            return;
        }
    }
}

bool OSC_Sim::checkFinished(){
    if(Error_found){
        cout << getId() << ": An error has been detected and the simulation will now end." << endl;
        return true;
    }
    if(Enable_exciton_diffusion_test){
        return (N_excitons_recombined==N_tests);
    }
    if(Enable_dynamics_test){
        if(getTime()>Dynamics_transient_end || (N_excitons==0 && N_electrons==0 && N_holes==0)){
            return true;
        }
        return false;
    }
    if(Enable_ToF_test){
        return (N_electrons_collected==N_tests || N_holes_collected==N_tests || N_electrons_created>2*N_tests || N_holes_created>2*N_tests);
    }
    if(Enable_IQE_test){
        if(isLightOn && N_excitons_created==N_tests){
            removeEvent(&exciton_creation_event);
            isLightOn = false;
        }
        if(N_excitons_created==N_tests && N_excitons==0 && N_electrons==0 && N_holes==0){
            return true;
        }
        if(getTime()>IQE_time_cutoff){
            return true;
        }
        return false;
    }
    cout << getId() << ": Error checking simulation finish conditions.  The simulation will now end." << endl;
    return true;
}

bool OSC_Sim::createImportedMorphology(){
    string file_info;
    string line;
    stringstream ss;
    Coords coords;
    short type = 0;
    int length,width,height;
    int site_count = 0;
    // Get input morphology file information
    getline(*Morphology_file,line);
    file_info = line;
    // Parse file info
    if(file_info.compare("Ising_OPV v3.2 - compressed format")!=0){
        cout << getId() << ": Error! Morphology file format not recognized. Only compressed morphologies created using Ising_OPV v3.2 are currently supported." << endl;
        Error_found = true;
        return false;
    }
    getline(*Morphology_file,line);
    length = atoi(line.c_str());
    getline(*Morphology_file,line);
    width = atoi(line.c_str());
    getline(*Morphology_file,line);
    height = atoi(line.c_str());
    if(lattice.getLength()!=length || lattice.getWidth()!=width || lattice.getHeight()!=height){
        cout << getId() << ": Error! Morphology lattice dimensions do not match the lattice dimensions defined in the parameter file." << endl;
        Error_found = true;
        return false;
    }
    // Skip 3 lines (domain size1, domain size2, blend ratio)
    getline(*Morphology_file,line);
    getline(*Morphology_file,line);
    getline(*Morphology_file,line);
    // Begin parsing morphology site data
    for(int x=0;x<lattice.getLength();x++){
        for(int y=0;y<lattice.getWidth();y++){
            for(int z=0;z<lattice.getHeight();z++){
                if(site_count==0){
                    if(!(*Morphology_file).good()){
                        cout << "Error parsing file.  End of file reached before expected." << endl;
                        Error_found = true;
                        return false;
                    }
                    getline(*Morphology_file,line);
                    type = (short)atoi(line.substr(0,1).c_str());
                    site_count = atoi(line.substr(1).c_str());
                }
                coords.setXYZ(x,y,z);
                sites[lattice.getSiteIndex(coords)].setType(type);
                if(type==(short)1){
                    N_donor_sites++;
                }
                else if(type==(short)2){
                    N_acceptor_sites++;
                }
                site_count--;
            }
        }
    }
    // Check for unassigned sites
    for(auto it=sites.begin();it!=sites.end();++it){
        if(it->getType()==(short)0){
            cout << getId() << ": Error! Unassigned site found after morphology import. Check the morphology file for errors." << endl;
            Error_found = true;
            return false;
        }
    }
    return true;
}

void OSC_Sim::deleteObject(const list<Object*>::iterator object_it){
    if((*object_it)->getName().compare(Exciton::name)==0){
        auto exciton_it = getExcitonIt(*object_it);
        // Remove the object from Simulation
        removeObject(object_it);
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
        auto polaron_it = getPolaronIt(*object_it);
        // Remove the object from Simulation
        removeObject(object_it);
        // Electron
        if(!(polaron_it->getCharge())){
            // Locate corresponding recombination event
            auto recombination_list_it = polaron_recombination_events.begin();
            advance(recombination_list_it,std::distance(electrons.begin(),polaron_it));
            // Locate corresponding hop event
            auto hop_list_it = electron_hop_events.begin();
            advance(hop_list_it,std::distance(electrons.begin(),polaron_it));
            // Locate corresponding dissociation event
            auto extraction_list_it = electron_extraction_events.begin();
            advance(extraction_list_it,std::distance(electrons.begin(),polaron_it));
            // Delete electron
            electrons.erase(polaron_it);
            // Delete polaron recombination event
            polaron_recombination_events.erase(recombination_list_it);
            // Delete electron hop event
            electron_hop_events.erase(hop_list_it);
            // Delete electron extraction event
            electron_extraction_events.erase(extraction_list_it);
        }
        // Hole
        else{
            // Locate corresponding hop event
            auto hop_list_it = hole_hop_events.begin();
            advance(hop_list_it,std::distance(holes.begin(),polaron_it));
            // Locate corresponding dissociation event
            auto extraction_list_it = hole_extraction_events.begin();
            advance(extraction_list_it,std::distance(holes.begin(),polaron_it));
            // Delete hole
            holes.erase(polaron_it);
            // Delete hole hop event
            hole_hop_events.erase(hop_list_it);
            // Delete hole extraction event
            hole_extraction_events.erase(extraction_list_it);
        }
    }
}

bool OSC_Sim::executeExcitonCreation(const list<Event*>::iterator event_it){
    // Determine coordinates for the new exciton
    const Coords coords_new = calculateExcitonCreationCoords();
    generateExciton(coords_new);
    // Find all nearby excitons and calculate their events
    auto neighbors = findRecalcNeighbors(coords_new);
    calculateObjectListEvents(neighbors);
    // Calculate next exciton creation event
    exciton_creation_event.calculateExecutionTime(R_exciton_generation_donor+R_exciton_generation_acceptor,getTime());
    return true;
}

bool OSC_Sim::executeExcitonDissociation(const list<Event*>::iterator event_it){
    // Get event info
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    Coords coords_dest = (*event_it)->getDestCoords();
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

bool OSC_Sim::executeExcitonHop(const list<Event*>::iterator event_it){
    if(lattice.isOccupied((*event_it)->getDestCoords())){
        cout << getId() << ": Error! Exciton hop cannot be executed. Destination site is already occupied." << endl;
        Error_found = true;
        return false;
    }
    else{
        if(loggingEnabled()){
            *Logfile << "Exciton " << (*((*event_it)->getObjectIt()))->getTag() << " hopping to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
        }
        return executeObjectHop(event_it);
    }
}

bool OSC_Sim::executeExcitonRecombine(const list<Event*>::iterator event_it){
    // Get event info
    int exciton_tag = (*((*event_it)->getObjectIt()))->getTag();
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
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
        *Logfile << "Exciton " << exciton_tag << " recombined at site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
    }
    // Find all nearby excitons and calculate their events
    auto neighbors = findRecalcNeighbors(coords_initial);
    calculateObjectListEvents(neighbors);
    return true;
}

bool OSC_Sim::executeNextEvent(){
    auto event_it = chooseNextEvent();
    if(*event_it==nullptr){
        cout << getId() << ": Error! The simulation has no events to execute." << endl;
        Error_found = true;
        return false;
    }
    string event_name = (*event_it)->getName();
    if(loggingEnabled()){
        *Logfile << "Executing " << event_name << " event" << endl;
    }
    // Update simulation time
    updateTime((*event_it)->getExecutionTime());
    // Update dynamics data
    if(Enable_dynamics_test){
        updateDynamicsData();
    }
    // Execute the chosen event
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
        cout << getId() << ": Error! Valid event not found when calling executeNextEvent" << endl;
        Error_found = true;
        return false;
    }
}

bool OSC_Sim::executeObjectHop(const list<Event*>::iterator event_it){
    // Get event info
    auto object_it = (*event_it)->getObjectIt();
    Coords coords_initial = (*object_it)->getCoords();
    Coords coords_dest = (*event_it)->getDestCoords();
    // Move the object in the Simulation
    moveObject((*event_it)->getObjectIt(),coords_dest);
    // Perform ToF test analysis
    if(Enable_ToF_test){
        updateToFData(object_it);
        // Check ToF time limit
        if((getTime()-(*object_it)->getCreationTime())>ToF_transient_end){
            // Update polaron counters
            if(!getObjectCharge(object_it)){
                N_electrons--;
            }
            else{
                N_holes--;
            }
            // Delete polaron and its events
            deleteObject(object_it);
        }
    }
    // Check if new ToF polarons need to be created
    if(Enable_ToF_test && N_holes==0 && N_electrons==0 && !checkFinished()){
        generateToFPolarons();
    }
    // Find all nearby objects and calculate their events
    else{
        auto neighbors = findRecalcNeighbors(coords_initial);
        auto neighbors2 = findRecalcNeighbors(coords_dest);
        neighbors.insert(neighbors.end(),neighbors2.begin(),neighbors2.end());
        removeObjectItDuplicates(neighbors);
        calculateObjectListEvents(neighbors);
    }
    return true;
}

bool OSC_Sim::executePolaronExtraction(const list<Event*>::iterator event_it){
    // Get event info
    bool charge = getObjectCharge((*event_it)->getObjectIt());
    int polaron_tag = (*((*event_it)->getObjectIt()))->getTag();
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
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
        if(!charge){
            *Logfile << "Electron " << polaron_tag << " was extracted from site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
        }
        else{
            *Logfile << "Hole " << polaron_tag << " was extracted from site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
        }
    }
    if(Enable_ToF_test && N_holes==0 && N_electrons==0 && !checkFinished()){
        generateToFPolarons();
    }
    else{
        // Find all nearby objects and calculate their events
        auto neighbors = findRecalcNeighbors(coords_initial);
        calculateObjectListEvents(neighbors);
    }
    return true;
}

bool OSC_Sim::executePolaronHop(const list<Event*>::iterator event_it){
    if(lattice.isOccupied((*event_it)->getDestCoords())){
        cout << getId() << ": Error! Polaron hop cannot be executed. Destination site is already occupied." << endl;
        Error_found = true;
        return false;
    }
    else{
        auto object_it = (*event_it)->getObjectIt();
        // Log event
        if(loggingEnabled()){
            if(!getObjectCharge(object_it)){
                *Logfile << "Electron " << (*object_it)->getTag() << " hopping to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
            }
            else{
                *Logfile << "Hole " << (*object_it)->getTag() << " hopping to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
            }
        }
        return executeObjectHop(event_it);
    }
}

bool OSC_Sim::executePolaronRecombination(const list<Event*>::iterator event_it){
    // Get event info
    auto object_it = (*event_it)->getObjectIt();
    int polaron_tag = (*object_it)->getTag();
    int target_tag = (*((*event_it)->getObjectTargetIt()))->getTag();
    Coords coords_initial = (*object_it)->getCoords();
    Coords coords_dest = (*event_it)->getDestCoords();
    // Delete polarons and their events
    deleteObject((*event_it)->getObjectTargetIt());
    deleteObject(object_it);
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
        *Logfile << "Electron " << polaron_tag << " at site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << " recombined with hole " << target_tag << " at site " << coords_dest.x << "," << coords_dest.y << "," << coords_dest.z << "." << endl;
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

void OSC_Sim::generateExciton(const Coords& coords){
    // Create the new exciton and add it to the simulation
    Exciton exciton_new(getTime(),N_excitons_created+1,coords);
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
    if(getSiteType(coords)==(short)1){
        N_excitons_created_donor++;
    }
    else{
        N_excitons_created_acceptor++;
    }
    N_excitons_created++;
    N_excitons++;
    // Log event
    if(loggingEnabled()){
        *Logfile << "Created exciton " << exciton_new.getTag() << " at site " << coords.x << "," << coords.y << "," << coords.z << "." << endl;
    }
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
        *Logfile << "Created electron " << electron_new.getTag() << " at site " << coords.x << "," << coords.y << "," << coords.z << "." << endl;
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
    N_holes_created++;
    N_holes++;
    // Log event
    if(loggingEnabled()){
        *Logfile << "Created hole " << hole_new.getTag() << " at site " << coords.x << "," << coords.y << "," << coords.z << "." << endl;
    }
}

void OSC_Sim::generateDynamicsExcitons(){
    int num = 0;
    Coords coords;
    int initial_excitons = ceil(Dynamics_initial_exciton_conc*intpow(1e-7*lattice.getUnitSize(),3)*lattice.getLength()*lattice.getWidth()*lattice.getHeight());
    cout << getId() << ": Generating " << initial_excitons << " initial excitons." << endl;
    while(num<initial_excitons){
        coords = calculateExcitonCreationCoords();
        generateExciton(coords);
        num++;
    }
    auto object_its = getAllObjectIts();
    calculateObjectListEvents(object_its);
}

void OSC_Sim::generateToFPolarons(){
    int num = 0;
    Coords coords;
    // Create electrons at the top plane of the lattice
    if(!ToF_polaron_type){
        coords.z = lattice.getHeight()-1;
    }
    // Create holes at the bottom plane of the lattice
    else{
        coords.z = 0;
    }
    ToF_start_times.assign(ToF_initial_polarons,getTime());
    ToF_start_positions.assign(ToF_initial_polarons,coords.z);
    ToF_start_energies.assign(ToF_initial_polarons,0);
    ToF_index_prev.assign(ToF_initial_polarons,0);
    auto energy_it = ToF_start_energies.begin();
    while(num<ToF_initial_polarons){
        coords.x = lattice.getRandomX();
        coords.y = lattice.getRandomY();
        // If the site is already occupied, pick a new site
        if(lattice.isOccupied(coords)){
            continue;
        }
        // If phase restriction is enabled, electrons cannot be created on donor sites
        else if(Enable_phase_restriction && !ToF_polaron_type && getSiteType(coords)==(short)1){
            continue;
        }
        // If phase restriction is enabled, holes cannot be created on acceptor sites
        else if(Enable_phase_restriction && ToF_polaron_type && getSiteType(coords)==(short)2){
            continue;
        }
        else if(!ToF_polaron_type){
            generateElectron(coords);
        }
        else{
            generateHole(coords);
        }
        *energy_it = getSiteEnergy(coords);
        energy_it++;
        num++;
    }
    auto object_its = getAllObjectIts();
    calculateObjectListEvents(object_its);
}

vector<double> OSC_Sim::getDiffusionData(){
    return diffusion_distances;
}

vector<int> OSC_Sim::getDynamicsTransientExcitons(){
    return transient_excitons;
}

vector<int> OSC_Sim::getDynamicsTransientElectrons(){
    return transient_electrons;
}

vector<int> OSC_Sim::getDynamicsTransientHoles(){
    return transient_holes;
}

vector<double> OSC_Sim::getDynamicsTransientTimes(){
    return transient_times;
}

list<Exciton>::iterator OSC_Sim::getExcitonIt(const Object* object_ptr){
    for(auto it=excitons.begin();it!=excitons.end();++it){
        if(object_ptr->getTag()==it->getTag()){
            return it;
        }
    }
}

int OSC_Sim::getN_bimolecular_recombinations(){
    return N_bimolecular_recombinations;
}

int OSC_Sim::getN_electrons_collected(){
    return N_electrons_collected;
}

int OSC_Sim::getN_electrons_created(){
    return N_electrons_created;
}

int OSC_Sim::getN_electrons_recombined(){
    return N_electrons_recombined;
}

int OSC_Sim::getN_excitons_created(){
    return N_excitons_created;
}

int OSC_Sim::getN_excitons_created(const short site_type){
    if(site_type==(short)1){
        return N_excitons_created_donor;
    }
    else{
        return N_excitons_created_acceptor;
    }
}

int OSC_Sim::getN_excitons_dissociated(){
    return N_excitons_dissociated;
}

int OSC_Sim::getN_excitons_recombined(){
    return N_excitons_recombined;
}

int OSC_Sim::getN_geminate_recombinations(){
    return N_geminate_recombinations;
}

int OSC_Sim::getN_holes_collected(){
    return N_holes_collected;
}

int OSC_Sim::getN_holes_created(){
    return N_holes_created;
}

int OSC_Sim::getN_holes_recombined(){
    return N_holes_recombined;
}

bool OSC_Sim::getObjectCharge(const list<Object*>::iterator object_it){
    if((*object_it)->getName().compare(Polaron::name)==0){
        return static_cast<Polaron*>(*object_it)->getCharge();
    }
}

list<Polaron>::iterator OSC_Sim::getPolaronIt(Object* object_ptr){
    if(object_ptr->getName().compare(Polaron::name)==0){
        // electrons
        if(!(static_cast<Polaron*>(object_ptr)->getCharge())){
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
    return sites[lattice.getSiteIndex(coords)].getEnergy();
}

short OSC_Sim::getSiteType(const Coords& coords){
    return static_cast<Site_OSC*>(*lattice.getSiteIt(coords))->getType();
}

vector<int> OSC_Sim::getToFTransientCounts(){
    return transient_counts;
}

vector<double> OSC_Sim::getToFTransientEnergies(){
    return transient_energies;
}

vector<double> OSC_Sim::getToFTransientTimes(){
    return transient_times;
}

vector<double> OSC_Sim::getToFTransientVelocities(){
    return transient_velocities;
}

vector<double> OSC_Sim::getTransitTimeData() {
	return transit_times;
}

void OSC_Sim::initializeArchitecture() {
	bool success;
	N_donor_sites = 0;
	N_acceptor_sites = 0;
	if (Enable_neat) {
		N_donor_sites = lattice.getNumSites();
		N_acceptor_sites = 0;
		for (auto it = sites.begin(); it != sites.end(); ++it) {
			it->setType(1);
		}
	}
	else if (Enable_bilayer) {
		Coords coords;
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					coords.setXYZ(x, y, z);
					if (z < Thickness_acceptor) {
						sites[lattice.getSiteIndex(coords)].setType(2);
						N_acceptor_sites++;
					}
					else {
						sites[lattice.getSiteIndex(coords)].setType(1);
						N_donor_sites++;
					}
				}
			}
		}
	}
	else if (Enable_random_blend) {
		vector<short> site_types;
		site_types.assign(lattice.getNumSites(), 1);
		for (int i = 0; i < (int)lattice.getNumSites()*Acceptor_conc; i++) {
			site_types[i] = 2;
			N_acceptor_sites++;
		}
		N_donor_sites = lattice.getNumSites() - N_acceptor_sites;
		shuffle(site_types.begin(), site_types.end(), gen);
		for (int i = 0; i < (int)sites.size(); i++) {
			sites[i].setType(site_types[i]);
		}
	}
	else if (Enable_import_morphology) {
		success = createImportedMorphology();
		if (!success) {
			return;
		}
	}
	if (Enable_gaussian_dos) {
		site_energies_donor.assign(N_donor_sites, 0);
		site_energies_acceptor.assign(N_acceptor_sites, 0);
		createGaussianDOSVector(site_energies_donor, 0, Energy_stdev_donor, gen);
		createGaussianDOSVector(site_energies_acceptor, 0, Energy_stdev_acceptor, gen);
	}
	else if (Enable_exponential_dos) {
		site_energies_donor.assign(N_donor_sites, 0);
		site_energies_acceptor.assign(N_acceptor_sites, 0);
		createExponentialDOSVector(site_energies_donor, 0, Energy_urbach_donor, gen);
		createExponentialDOSVector(site_energies_acceptor, 0, Energy_urbach_acceptor, gen);
	}
	else {
		site_energies_donor.push_back(0);
		site_energies_acceptor.push_back(0);
	}
	int donor_count = 0;
	int acceptor_count = 0;
	for (int i = 0; i < lattice.getNumSites(); i++) {
		if (Enable_gaussian_dos || Enable_exponential_dos) {
			if (sites[i].getType() == (short)1) {
				sites[i].setEnergyIt(site_energies_donor.begin() + donor_count);
				donor_count++;
			}
			else if (sites[i].getType() == (short)2) {
				sites[i].setEnergyIt(site_energies_acceptor.begin() + acceptor_count);
				acceptor_count++;
			}
			else {
				cout << getId() << ": Error! Undefined site type detected while assigning site energies." << endl;
				Error_found = true;
			}
		}
		else {
			if (sites[i].getType() == (short)1) {
				sites[i].setEnergyIt(site_energies_donor.begin());
			}
			else if (sites[i].getType() == (short)2) {
				sites[i].setEnergyIt(site_energies_acceptor.begin());
			}
			else {
				cout << getId() << ": Error! Undefined site type detected while assigning site energies." << endl;
				Error_found = true;
			}
		}
	}
	vector<Site*> site_ptrs((int)sites.size());
	for (int i = 0; i < (int)sites.size(); i++){
		site_ptrs[i] = &sites[i];
	}
	lattice.setSitePointers(site_ptrs);
}

void OSC_Sim::outputStatus(){
    cout << getId() << ": Time = " << getTime() << " seconds.\n";
    if(Enable_ToF_test){
        if(!ToF_polaron_type){
            cout << getId() << ": " << N_electrons_collected << " out of " << N_electrons_created << " electrons have been collected and " << getN_events_executed() << " events have been executed.\n";
            cout << getId() << ": There are currently " << N_electrons << " electrons in the lattice:\n";
            for(auto it=electrons.begin();it!=electrons.end();++it){
                cout << getId() << ": Electron " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
            }
        }
        else{
            cout << getId() << ": " << N_holes_collected << " out of " << N_holes_created << " holes have been collected and " << getN_events_executed() << " events have been executed.\n";
            cout << getId() << ": There are currently " << N_holes << " holes in the lattice:\n";
            for(auto it=holes.begin();it!=holes.end();++it){
                cout << getId() << ": Hole " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
            }
        }
    }
    if(Enable_exciton_diffusion_test){
        cout << getId() << ": " << N_excitons_created << " excitons have been created and " << getN_events_executed() << " events have been executed.\n";
        cout << getId() << ": There are currently " << N_excitons << " excitons in the lattice:\n";
        for(auto it=excitons.begin();it!=excitons.end();++it){
            cout << getId() << ": Exciton " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
        }
    }
    if(Enable_IQE_test || Enable_dynamics_test){
        cout << getId() << ": " << N_excitons_created << " excitons have been created and " << getN_events_executed() << " events have been executed.\n";
        cout << getId() << ": There are currently " << N_excitons << " excitons in the lattice:\n";
        for(auto it=excitons.begin();it!=excitons.end();++it){
            cout << getId() << ": Exciton " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
        }
        cout << getId() << ": There are currently " << N_electrons << " electrons in the lattice:\n";
        for(auto it=electrons.begin();it!=electrons.end();++it){
            cout << getId() << ": Electron " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
        }
        cout << getId() << ": There are currently " << N_holes << " holes in the lattice:\n";
        for(auto it=holes.begin();it!=holes.end();++it){
            cout << getId() << ": Hole " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
        }
    }
    cout.flush();
}

bool OSC_Sim::siteContainsHole(const Coords& coords){
    auto object_it = (*lattice.getSiteIt(coords))->getObjectIt();
    if((*object_it)->getName().compare(Polaron::name)==0){
        return static_cast<Polaron*>(*object_it)->getCharge();
    }
    return false;
}

void OSC_Sim::updateDynamicsData(){
    static const double step_size = 1.0/(double)Dynamics_pnts_per_decade;
    static int index_prev = -1;
    // Only start dynamics transient calculations if the time elapsed since exciton creation is larger than the transient start time
    if(getTime()>Dynamics_transient_start){
        // If enough time has passed, output next timestep data
        if(getTime()>transient_times[index_prev+1]){
            int index = floor((log10(getTime())-log10(Dynamics_transient_start))/step_size);
            while(index!=0 && index_prev<index-1 && index_prev+1<(int)transient_times.size()){
                transient_excitons[index_prev+1] = N_excitons;
                transient_electrons[index_prev+1] = N_electrons;
                transient_holes[index_prev+1] = N_holes;
                index_prev++;
            }
            if(index<(int)transient_times.size()){
                transient_excitons[index] = N_excitons;
                transient_electrons[index] = N_electrons;
                transient_holes[index] = N_holes;
                index_prev = index;
            }
        }
    }
}

void OSC_Sim::updateToFData(const list<Object*>::iterator object_it){
    static const double step_size = 1.0/(double)ToF_pnts_per_decade;
    // Only start ToF transient calculations if the time elapsed since polaron creation is larger than the ToF transient start time
    if((getTime()-(*object_it)->getCreationTime())>ToF_transient_start){
        // Get polaron and previous timestep info
        auto polaron_it = getPolaronIt(*object_it);
        auto start_time_it = ToF_start_times.begin();
        if(!polaron_it->getCharge()){
            advance(start_time_it,distance(electrons.begin(),polaron_it));
        }
        else{
            advance(start_time_it,distance(holes.begin(),polaron_it));
        }
        // If enough time has passed, output next timestep data
        if(log10(getTime()-(*object_it)->getCreationTime())-log10(*start_time_it-(*object_it)->getCreationTime())>step_size){
            int index = floor((log10(getTime()-(*object_it)->getCreationTime())-log10(ToF_transient_start))/step_size);
            // Get polaron site energy for previous timestep
            auto start_energy_it = ToF_start_energies.begin();
            auto index_prev_it = ToF_index_prev.begin();
            if(!polaron_it->getCharge()){
                advance(start_energy_it,distance(electrons.begin(),polaron_it));
                advance(index_prev_it,distance(electrons.begin(),polaron_it));
            }
            else{
                advance(start_energy_it,distance(holes.begin(),polaron_it));
                advance(index_prev_it,distance(holes.begin(),polaron_it));
            }
            int index_prev = *index_prev_it;
            while(index!=0 && index_prev<index-1 && index_prev+1<(int)transient_times.size()){
                transient_counts[index_prev+1]++;
                // transient_velocities[index_prev+1] += 0;
                transient_energies[index_prev+1] += *start_energy_it;
                *start_time_it += pow(10,log10(*start_time_it)+step_size);
                index_prev++;
            }
            // Get polaron position for previous timestep
            auto start_position_it = ToF_start_positions.begin();
            if(!polaron_it->getCharge()){
                advance(start_position_it,distance(electrons.begin(),polaron_it));
            }
            else{
                advance(start_position_it,distance(holes.begin(),polaron_it));
            }
            if(index<(int)transient_times.size()){
                transient_counts[index]++;
                transient_velocities[index] += (1e-7*lattice.getUnitSize()*abs((*object_it)->getCoords().z - *start_position_it))/(getTime()- *start_time_it);
                transient_energies[index] += getSiteEnergy((*object_it)->getCoords());
                *start_time_it = getTime();
                *start_position_it = (*object_it)->getCoords().z;
                *start_energy_it = getSiteEnergy((*object_it)->getCoords());
                *index_prev_it = index;
            }
        }

    }
}
