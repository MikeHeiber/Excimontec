// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef OSC_SIM_H
#define OSC_SIM_H

#include "KMC_Lattice/Utils.h"
#include "KMC_Lattice/Simulation.h"
#include "KMC_Lattice/Object.h"
#include "KMC_Lattice/Event.h"
#include "Exciton.h"
#include "Polaron.h"
#include <algorithm>

using namespace std;

struct Parameters_OPV : Parameters_Simulation{
    // Additional General Parameters
    double Bias;
    // Morphology Parameters
    bool Enable_neat; // Neat takes on donor properties
    bool Enable_bilayer;
    int Thickness_donor; // sites
    int Thickness_acceptor; // sites
    bool Enable_random_blend;
    double Acceptor_conc;
    // Test Parameters
    int N_tests;
    bool Enable_exciton_diffusion_test;
    bool Enable_ToF_test;
    bool ToF_polaron_type;
    int ToF_initial_polarons;
    double ToF_transient_start;
    double ToF_transient_end;
    int ToF_pnts_per_decade;
    bool Enable_IQE_test;
    // Exciton Parameters
    double Exciton_generation_rate_donor;
    double Exciton_generation_rate_acceptor;
    double Exciton_lifetime_donor; // seconds
    double Exciton_lifetime_acceptor; // seconds
    double R_exciton_hopping_donor;
    double R_exciton_hopping_acceptor;
    int FRET_cutoff;
    double E_exciton_binding_donor;
    double E_exciton_binding_acceptor;
    double R_exciton_dissociation_donor;
    double R_exciton_dissociation_acceptor;
    int Exciton_dissociation_cutoff; // nm
    // Polaron Parameters
    bool Enable_phase_restriction;
    double R_polaron_hopping_donor;
    double R_polaron_hopping_acceptor;
    double Polaron_localization_donor; // nm^-1
    double Polaron_localization_acceptor; // nm^-1
    bool Enable_miller_abrahams;
    bool Enable_marcus;
    double Reorganization_donor;
    double Reorganization_acceptor;
    double R_polaron_recombination;
    int Polaron_hopping_cutoff; // nm
    // Additional Lattice Parameters
    double Homo_donor;
    double Lumo_donor;
    double Homo_acceptor;
    double Lumo_acceptor;
    bool Enable_gaussian_dos;
    double Energy_stdev_donor; // eV
    double Energy_stdev_acceptor; // eV
    bool Enable_exponential_dos;
    double Energy_urbach_donor; // eV
    double Energy_urbach_acceptor; // eV
    // Coulomb Calculation Parameters
    double Dielectric_donor;
    double Dielectric_acceptor;
    int Coulomb_cutoff; // nm
};

class Site_OSC : public Site{
    public:
        double getEnergy(){return *energy_it;}
        short getType(){return type;}
        void setEnergyIt(const vector<double>::iterator it){energy_it = it;}
        void setType(const short site_type){type = site_type;}
    private:
        vector<double>::iterator energy_it;
        short type; //  type 1 represent donor, type 2 represents acceptor
};

class OSC_Sim : public Simulation{
    public:
        // Functions
        OSC_Sim(const Parameters_OPV& params,const int id);
        double calculateDiffusionLength_avg();
        double calculateDiffusionLength_stdev();
        vector<double> calculateTransitTimeDist(const vector<double>& data,const int counts);
        double calculateTransitTime_avg();
        double calculateTransitTime_stdev();
        double calculateMobility_avg();
        double calculateMobility_stdev();
        bool checkFinished();
        bool executeNextEvent();
        vector<double> getDiffusionData();
        vector<int> getToFTransientCounts();
        vector<double> getToFTransientEnergies();
        vector<double> getToFTransientTimes();
        vector<double> getToFTransientVelocities();
        vector<double> getTransitTimeData();
        int getN_excitons_created();
        int getN_excitons_dissociated();
        int getN_excitons_recombined();
        int getN_electrons_created();
        int getN_electrons_collected();
        int getN_electrons_recombined();
        int getN_holes_created();
        int getN_holes_collected();
        int getN_holes_recombined();
        int getN_geminate_recombinations();
        int getN_bimolecular_recombinations();
        void outputStatus();
    protected:

    private:
        // Additional General Parameters
        double Bias;
        // Morphology Parameters
        bool Enable_neat; // Neat takes on donor properties
        bool Enable_bilayer;
        int Thickness_donor; // sites
        int Thickness_acceptor; // sites
        bool Enable_random_blend;
        double Acceptor_conc;
        // Test Parameters
        int N_tests;
        bool Enable_exciton_diffusion_test;
        bool Enable_ToF_test;
        bool ToF_polaron_type;
        int ToF_initial_polarons;
        double ToF_transient_start;
        double ToF_transient_end;
        int ToF_pnts_per_decade;
        bool Enable_IQE_test;
        // Exciton Parameters
        double Exciton_generation_rate_donor;
        double Exciton_generation_rate_acceptor;
        double Exciton_lifetime_donor; // seconds
        double Exciton_lifetime_acceptor; // seconds
        double R_exciton_hopping_donor;
        double R_exciton_hopping_acceptor;
        int FRET_cutoff;
        double E_exciton_binding_donor;
        double E_exciton_binding_acceptor;
        double R_exciton_dissociation_donor;
        double R_exciton_dissociation_acceptor;
        int Exciton_dissociation_cutoff; // nm
        // Polaron Parameters
        bool Enable_phase_restriction;
        double R_polaron_hopping_donor;
        double R_polaron_hopping_acceptor;
        double Polaron_localization_donor; // nm^-1
        double Polaron_localization_acceptor; // nm^-1
        bool Enable_miller_abrahams;
        bool Enable_marcus;
        double Reorganization_donor;
        double Reorganization_acceptor;
        double R_polaron_recombination;
        int Polaron_hopping_cutoff; // nm
        // Additional Lattice Parameters
        double Homo_donor;
        double Lumo_donor;
        double Homo_acceptor;
        double Lumo_acceptor;
        bool Enable_gaussian_dos;
        double Energy_stdev_donor; // eV
        double Energy_stdev_acceptor; // eV
        bool Enable_exponential_dos;
        double Energy_urbach_donor;
        double Energy_urbach_acceptor;
        // Coulomb Calculation Parameters
        double Dielectric_donor;
        double Dielectric_acceptor;
        int Coulomb_cutoff; // nm
        // Additional Output Files
        //
        // Additional Parameters
        bool Error_found;
        double R_exciton_generation_donor;
        double R_exciton_generation_acceptor;
        // Site Data Structure
        vector<Site_OSC> sites;
        // Object Data Structures
        list<Exciton> excitons;
        list<Polaron> electrons;
        list<Polaron> holes;
        // Event Data Structures
        list<Exciton_Hop> exciton_hop_events;
        list<Exciton_Recombination> exciton_recombination_events;
        list<Exciton_Dissociation> exciton_dissociation_events;
        list<Polaron_Hop> electron_hop_events;
        list<Polaron_Hop> hole_hop_events;
        list<Polaron_Recombination> polaron_recombination_events;
        list<Polaron_Extraction> electron_extraction_events;
        list<Polaron_Extraction> hole_extraction_events;
        // Additional Data Structures
        vector<double> Coulomb_table;
        vector<double> E_potential;
        vector<double> site_energies_donor;
        vector<double> site_energies_acceptor;
        vector<double> diffusion_distances;
        list<int> ToF_start_positions;
        list<int> ToF_index_prev;
        list<double> ToF_start_times;
        list<double> ToF_start_energies;
        vector<double> transient_times;
        vector<int> transient_counts;
        vector<double> transient_velocities;
        vector<double> transient_energies;
        vector<double> transit_times;
        Exciton_Creation exciton_creation_event;
        list<Event*>::iterator exciton_creation_it;
        // Additional Counters
        int N_donor_sites;
        int N_acceptor_sites;
        int N_excitons_created;
        int N_excitons_recombined;
        int N_excitons_dissociated;
        int N_excitons_quenched;
        int N_excitons;
        int N_electrons_created;
        int N_electrons_recombined;
        int N_electrons_collected;
        int N_electrons;
        int N_holes_created;
        int N_holes_recombined;
        int N_holes_collected;
        int N_holes;
        int N_geminate_recombinations;
        int N_bimolecular_recombinations;
        int N_electron_surface_recombinations;
        int N_hole_surface_recombinations;
        // Additional Functions
        double calculateCoulomb(const list<Object*>::iterator object_it,const Coords& coords);
        double calculateCoulomb(const bool charge,const Coords& coords);
        Coords calculateExcitonCreationCoords();
        void calculateExcitonEvents(const list<Object*>::iterator object_it);
        void calculateObjectListEvents(const vector<list<Object*>::iterator>& object_it_vec);
        void calculatePolaronEvents(const list<Object*>::iterator object_it);
        void deleteObject(const list<Object*>::iterator object_it);
        // Exciton Event Execution Functions
        bool executeExcitonCreation(const list<Event*>::iterator event_it);
        bool executeExcitonHop(const list<Event*>::iterator event_it);
        bool executeExcitonRecombine(const list<Event*>::iterator event_it);
        bool executeExcitonDissociation(const list<Event*>::iterator event_it);
        bool executeExcitonIntersystemCrossing(const list<Event*>::iterator event_it);
        bool executeExcitonExcitonAnnihilation(const list<Event*>::iterator event_it);
        bool executeExcitonPolaronAnnihilation(const list<Event*>::iterator event_it);
        // General Event Functions
        bool executeObjectHop(const list<Event*>::iterator event_it);
        // Polaron Event Execution Functions
        bool executePolaronHop(const list<Event*>::iterator event_it);
        bool executePolaronRecombination(const list<Event*>::iterator event_it);
        bool executePolaronExtraction(const list<Event*>::iterator event_it);
        void generateElectron(const Coords& coords,int tag);
        void generateHole(const Coords& coords,int tag);
        void generateToFPolarons();
        list<Exciton>::iterator getExcitonIt(const Object* object_ptr);
        bool getObjectCharge(const list<Object*>::iterator object_it);
        list<Polaron>::iterator getPolaronIt(Object* object_ptr);
        double getSiteEnergy(const Coords& coords);
        short getSiteType(const Coords& coords);
        void initializeArchitecture();
        bool siteContainsHole(const Coords& coords);
        void updateToFData(const list<Object*>::iterator object_it);
};

#endif //OSC_SIM_H
