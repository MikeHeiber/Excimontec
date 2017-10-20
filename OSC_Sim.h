// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef OSC_SIM_H
#define OSC_SIM_H

#include "KMC_Lattice/Simulation.h"
#include "KMC_Lattice/Site.h"
#include "Exciton.h"
#include "Polaron.h"
#include <algorithm>
#include <numeric>

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
    bool Enable_import_morphology;
    std::ifstream* Morphology_file;
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
    double IQE_time_cutoff;
    bool Enable_dynamics_test;
    bool Enable_dynamics_extraction;
    double Dynamics_initial_exciton_conc;
    double Dynamics_transient_start;
    double Dynamics_transient_end;
    int Dynamics_pnts_per_decade;
    // Exciton Parameters
    double Exciton_generation_rate_donor;
    double Exciton_generation_rate_acceptor;
    double Singlet_lifetime_donor; // seconds
    double Singlet_lifetime_acceptor; // seconds
	double Triplet_lifetime_donor; // seconds
	double Triplet_lifetime_acceptor; // seconds
    double R_singlet_hopping_donor;
    double R_singlet_hopping_acceptor;
	double R_triplet_hopping_donor;
	double R_triplet_hopping_acceptor;
	double Triplet_localization_donor;
	double Triplet_localization_acceptor;
	bool Enable_FRET_triplet_annihilation;
	double R_exciton_exciton_annihilation_donor;
	double R_exciton_exciton_annihilation_acceptor;
	double R_exciton_polaron_annihilation_donor;
	double R_exciton_polaron_annihilation_acceptor;
    int FRET_cutoff;
    double E_exciton_binding_donor;
    double E_exciton_binding_acceptor;
    double R_exciton_dissociation_donor;
    double R_exciton_dissociation_acceptor;
    int Exciton_dissociation_cutoff; // nm
	double R_exciton_isc_donor;
	double R_exciton_isc_acceptor;
	double R_exciton_risc_donor;
	double R_exciton_risc_acceptor;
	double E_exciton_ST_donor;
	double E_exciton_ST_acceptor;
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
    bool Enable_gaussian_polaron_delocalization;
    double Polaron_delocalization_length;
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
	bool Enable_correlated_disorder;
	double Disorder_correlation_length; // nm
	bool Enable_gaussian_kernel;
	bool Enable_power_kernel;
	int Power_kernel_exponent; // must be negative
    // Coulomb Calculation Parameters
    double Dielectric_donor;
    double Dielectric_acceptor;
    int Coulomb_cutoff; // nm
};

class Site_OSC : public Site{
    public:
        double getEnergy() const{return *energy_it;}
        short getType() const{return type;}
		void setEnergy(const double energy) { *energy_it = energy; }
        void setEnergyIt(const std::vector<double>::iterator it){energy_it = it;}
        void setType(const short site_type){type = site_type;}
    private:
		std::vector<double>::iterator energy_it;
        short type = 0; //  type 1 represent donor, type 2 represents acceptor
};

class OSC_Sim : public Simulation{
    public:
        // Functions
		OSC_Sim();
		virtual ~OSC_Sim();
        bool init(const Parameters_OPV& params,const int id);
        double calculateDiffusionLength_avg() const;
        double calculateDiffusionLength_stdev() const;
		std::vector<std::pair<double,double>> calculateDOSCorrelation(const double cutoff_radius);
		std::vector<double> calculateTransitTimeDist(const std::vector<double>& data,const int counts) const;
        double calculateTransitTime_avg() const;
        double calculateTransitTime_stdev() const;
		std::vector<double> calculateMobilities(const std::vector<double>& transit_times) const;
        double calculateMobility_avg() const;
        double calculateMobility_stdev() const;
        bool checkFinished() const;
		bool checkParameters(const Parameters_OPV& params) const;
        bool executeNextEvent();
		std::vector<double> getDiffusionData() const;
		std::vector<std::pair<double, double>> getDOSCorrelationData() const;
		std::vector<int> getDynamicsTransientExcitons() const;
		std::vector<int> getDynamicsTransientElectrons() const;
		std::vector<int> getDynamicsTransientHoles() const;
		std::vector<double> getDynamicsTransientTimes() const;
		std::vector<double> getSiteEnergies(const short site_type) const;
		std::vector<std::string> getChargeExtractionMap(const bool charge) const;
		std::vector<int> getToFTransientCounts() const;
		std::vector<double> getToFTransientEnergies() const;
		std::vector<double> getToFTransientTimes() const;
		std::vector<double> getToFTransientVelocities() const;
		std::vector<double> getTransitTimeData() const;
        int getN_excitons_created() const;
        int getN_excitons_created(const short site_type) const;
        int getN_excitons_dissociated() const;
        int getN_excitons_recombined() const;
        int getN_electrons_created() const;
        int getN_electrons_collected() const;
        int getN_electrons_recombined() const;
        int getN_holes_created() const;
        int getN_holes_collected() const;
        int getN_holes_recombined() const;
        int getN_geminate_recombinations() const;
        int getN_bimolecular_recombinations() const;
		int getN_transient_cycles() const;
        void outputStatus();
		void reassignSiteEnergies();
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
        bool Enable_import_morphology;
		std::ifstream* Morphology_file;
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
        double IQE_time_cutoff;
        bool Enable_dynamics_test;
        bool Enable_dynamics_extraction;
        double Dynamics_initial_exciton_conc;
        double Dynamics_transient_start;
        double Dynamics_transient_end;
        int Dynamics_pnts_per_decade;
        // Exciton Parameters
        double Exciton_generation_rate_donor;
        double Exciton_generation_rate_acceptor;
        double Singlet_lifetime_donor; // seconds
        double Singlet_lifetime_acceptor; // seconds
		double Triplet_lifetime_donor; // seconds
		double Triplet_lifetime_acceptor; // seconds
        double R_singlet_hopping_donor;
        double R_singlet_hopping_acceptor;
		double R_triplet_hopping_donor;
		double R_triplet_hopping_acceptor;
		double Triplet_localization_donor;
		double Triplet_localization_acceptor;
		bool Enable_FRET_triplet_annihilation;
		double R_exciton_exciton_annihilation_donor;
		double R_exciton_exciton_annihilation_acceptor;
		double R_exciton_polaron_annihilation_donor;
		double R_exciton_polaron_annihilation_acceptor;
        int FRET_cutoff;
        double E_exciton_binding_donor;
        double E_exciton_binding_acceptor;
        double R_exciton_dissociation_donor;
        double R_exciton_dissociation_acceptor;
        int Exciton_dissociation_cutoff; // nm
		double R_exciton_isc_donor;
		double R_exciton_isc_acceptor;
		double R_exciton_risc_donor;
		double R_exciton_risc_acceptor;
		double E_exciton_ST_donor;
		double E_exciton_ST_acceptor;
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
        bool Enable_gaussian_polaron_delocalization;
        double Polaron_delocalization_length;
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
		bool Enable_correlated_disorder;
		double Disorder_correlation_length; // nm
		bool Enable_gaussian_kernel;
		bool Enable_power_kernel;
		int Power_kernel_exponent; // must be negative
        // Coulomb Calculation Parameters
        double Dielectric_donor;
        double Dielectric_acceptor;
        int Coulomb_cutoff; // nm
        // Additional Output Files
        //
        // Additional Parameters
        bool isLightOn;
        double R_exciton_generation_donor;
        double R_exciton_generation_acceptor;
		double ToF_creation_time;
		int ToF_index_prev;
        // Site Data Structure
		std::vector<Site_OSC> sites;
        // Object Data Structures
		std::list<Exciton> excitons;
		std::list<Polaron> electrons;
		std::list<Polaron> holes;
        // Event Data Structures
		std::list<Exciton_Hop> exciton_hop_events;
		std::list<Exciton_Recombination> exciton_recombination_events;
		std::list<Exciton_Dissociation> exciton_dissociation_events;
		std::list<Exciton_Exciton_Annihilation> exciton_exciton_annihilation_events;
		std::list<Exciton_Polaron_Annihilation> exciton_polaron_annihilation_events;
		std::list<Exciton_Intersystem_Crossing> exciton_intersystem_crossing_events;
		std::list<Polaron_Hop> electron_hop_events;
		std::list<Polaron_Hop> hole_hop_events;
		std::list<Polaron_Recombination> polaron_recombination_events;
		std::list<Polaron_Extraction> electron_extraction_events;
		std::list<Polaron_Extraction> hole_extraction_events;
        // Additional Data Structures
		std::vector<double> Coulomb_table;
		std::vector<double> E_potential;
		std::vector<double> site_energies_donor;
		std::vector<double> site_energies_acceptor;
		std::vector<std::pair<double, double>> DOS_correlation_data;
		std::vector<double> diffusion_distances;
		std::vector<int> ToF_polaron_tags;
		std::vector<int> ToF_positions_prev;
		std::vector<double> ToF_times_prev;
		std::vector<double> ToF_energies_prev;
		std::vector<int> Electron_extraction_data;
		std::vector<int> Hole_extraction_data;
		std::vector<double> transient_times;
		std::vector<int> transient_counts;
		std::vector<double> transient_velocities;
		std::vector<double> transient_energies;
		std::vector<double> transit_times;
		std::vector<int> transient_excitons;
		std::vector<int> transient_electrons;
		std::vector<int> transient_holes;
        Exciton_Creation exciton_creation_event;
		std::list<Event*>::iterator exciton_creation_it;
        // Additional Counters
        int N_donor_sites;
        int N_acceptor_sites;
        int N_excitons_created = 0;
        int N_excitons_created_donor = 0;
        int N_excitons_created_acceptor = 0;
        int N_excitons_recombined = 0;
        int N_excitons_dissociated = 0;
		int N_exciton_exciton_annihilations = 0;
		int N_exciton_polaron_annihilations = 0;
		int N_exciton_intersystem_crossings = 0;
		int N_exciton_reverse_intersystem_crossings = 0;
        int N_excitons_quenched = 0;
        int N_excitons = 0;
        int N_electrons_created = 0;
        int N_electrons_recombined = 0;
        int N_electrons_collected = 0;
        int N_electrons = 0;
        int N_holes_created = 0;
        int N_holes_recombined = 0;
        int N_holes_collected = 0;
        int N_holes = 0;
        int N_geminate_recombinations = 0;
        int N_bimolecular_recombinations = 0;
        int N_electron_surface_recombinations = 0;
        int N_hole_surface_recombinations = 0;
		int N_transient_cycles = 0;
        // Additional Functions
        double calculateCoulomb(const std::list<Polaron>::iterator,const Coords& coords) const;
        double calculateCoulomb(const bool charge,const Coords& coords) const;
        Coords calculateExcitonCreationCoords();
        void calculateExcitonEvents(Object* object_ptr);
        void calculateObjectListEvents(const std::vector<Object*>& object_ptr_vec);
        void calculatePolaronEvents(Object* object_ptr);
		void createCorrelatedDOS(const double correlation_length);
        bool createImportedMorphology();
        void deleteObject(Object* object_ptr);
        // Exciton Event Execution Functions
        bool executeExcitonCreation();
        bool executeExcitonHop(const std::list<Event*>::iterator event_it);
        bool executeExcitonRecombine(const std::list<Event*>::iterator event_it);
        bool executeExcitonDissociation(const std::list<Event*>::iterator event_it);
        bool executeExcitonIntersystemCrossing(const std::list<Event*>::iterator event_it);
        bool executeExcitonExcitonAnnihilation(const std::list<Event*>::iterator event_it);
        bool executeExcitonPolaronAnnihilation(const std::list<Event*>::iterator event_it);
        // General Event Functions
        bool executeObjectHop(const std::list<Event*>::iterator event_it);
        // Polaron Event Execution Functions
        bool executePolaronHop(const std::list<Event*>::iterator event_it);
        bool executePolaronRecombination(const std::list<Event*>::iterator event_it);
        bool executePolaronExtraction(const std::list<Event*>::iterator event_it);
        void generateExciton(const Coords& coords);
        void generateElectron(const Coords& coords,int tag);
        void generateHole(const Coords& coords,int tag);
        void generateDynamicsExcitons();
        void generateToFPolarons();
		std::list<Exciton>::iterator getExcitonIt(const Object* object_ptr);
		std::list<Polaron>::iterator getPolaronIt(const Object* object_ptr);
        double getSiteEnergy(const Coords& coords) const;
        short getSiteType(const Coords& coords) const;
        bool initializeArchitecture();
        bool siteContainsHole(const Coords& coords);
        void updateDynamicsData();
        void updateToFData();
};

#endif //OSC_SIM_H
