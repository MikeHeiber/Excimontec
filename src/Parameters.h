// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCIMONTEC_PARAMETERS_H
#define EXCIMONTEC_PARAMETERS_H

#include "Lattice.h"
#include "Parameters_Simulation.h"
#include "Utils.h"

namespace Excimontec {

	//! \brief This class extends the KMC_Lattice::Parameters_Simulation class to create a parameters object that holds all input parameters needed by a simulation.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2017-2019
	class Parameters : public KMC_Lattice::Parameters_Simulation {
	public:

		// main Parameters ---------------------------------------------------------------------------------------
		// These parameters are only used within main.cpp to setup the simulation.

		//! Specifies whether or not to import a single morphology file and use it on all processors.
		bool Enable_import_morphology_single = false;

		//! Specifies whether or not to import a set of morphologies to be assigned to different processors.
		bool Enable_import_morphology_set = false;

		//! \brief Defines the morphology set naming format.
		//! The standard naming format is morphology_#.txt where # is a wildcard for different numbers.
		std::string Morphology_set_format = "";

		//! Defines the number of morphologies to be tested from the set.
		int N_test_morphologies = 0;

		//! \brief Defines the number of morphologies in the set.
		//! It is assumed that morphologies in the set are number from 0 to N_morphology_set_size-1.
		int N_morphology_set_size = 0;

		//! Specifies whether or not extraction map data should be output to text file at the end of the simulation.
		bool Enable_extraction_map_output = false;

		// Additional General Parameters -------------------------------------------------------------------------

		//! Defines the internal electrical potential across the semiconductor layer
		double Internal_potential;

		// Morphology Parameters ---------------------------------------------------------------------------------

		//! \brief Specifies whether or not to use a neat single phase material device architecture.
		//! A neat material morphology uses the donor material properties.
		bool Enable_neat;

		//! Specifies whether or not to use a bilayer device architecture.
		bool Enable_bilayer;

		//! Defines the donor layer thickness in number of sites when using bilayer architecture.
		int Thickness_donor; // sites

		//! Defines the acceptor layer thickness in number of sites when using bilayer architecture.
		int Thickness_acceptor; // sites

		//! Specifies whether or not to use a random blend material device architecture.
		bool Enable_random_blend;

		//! Defines the acceptor site volume concentration when using the random blend device architecture.
		double Acceptor_conc;

		//! Specifies whether or not an imported morphology is being used to define the device architecture.
		bool Enable_import_morphology;

		//! The name of the morphology file to be imported.
		std::string Morphology_filename;

		// Simulation Test Parameters ----------------------------------------------------------------------------------
		// These parameters determine which simulation test will be run and the settings for the simulation test.

		//! Used as a generic counter for determining when a given simulation test is completed.
		int N_tests;

		//! Specifies whether or not the exciton diffusion test is to be run.
		bool Enable_exciton_diffusion_test;

		//! Specifies whether or not to run the time-of-flight charge transport test.
		bool Enable_ToF_test;

		//! \brief Specifies the type of polaron to test when performing the time-of-flight charge transport test. 
		//! (false for electron and true for hole)
		bool ToF_polaron_type;

		//! Defines the number of initial polarons to place in the simulation for the time-of-flight charge transport test.
		int ToF_initial_polarons;

		//! Specifies whether or not to use the random placement method for creating the initial polarons in the time-of-flight charge transport test.
		bool Enable_ToF_random_placement;

		//! Specifies whether or not to use the energy placement method for creating the initial polarons in the time-of-flight charge transport test.
		bool Enable_ToF_energy_placement;

		//! Defines the preferred energy within the DOS for initial polaron placement when using the energy placement method with the time-of-flight charge transport test.
		double ToF_placement_energy;

		//! Defines the start time of the transient data collected during the time-of-flight charge transport test.
		double ToF_transient_start;

		//! Defines the end time of the transient data collected during the time-of-flight charge transport test.
		double ToF_transient_end;

		//! Defines the collected transient data time resolution in number of points per decade for the time-of-flight charge transport test.
		int ToF_pnts_per_decade;

		//! Specifies whether or not to run the internal quantum efficiency test.
		bool Enable_IQE_test;

		//! Defines the simulation time cutoff for the internal quantum efficiency test.
		double IQE_time_cutoff;

		//! Specifies whether or not to run the dynamics test.
		bool Enable_dynamics_test;

		//! Specifies whether or not to allow polaron extraction at the electrodes during the dynamics test.
		bool Enable_dynamics_extraction;

		//! Defines the initial exciton concentration used by the dynamics test.
		double Dynamics_initial_exciton_conc;

		//! Defines the start time of the transient data collected during the dynamics test.
		double Dynamics_transient_start;

		//! Defines the end time of the transient data collected during the dynamics test.
		double Dynamics_transient_end;

		//! Defines the collected transient data time resolution in number of points per decade for the dynamics test.
		int Dynamics_pnts_per_decade;

		//! Specifies whether or not the run the steady state charge transport test.
		bool Enable_steady_transport_test;

		//! Defines the steady state charge carrier density used by the steady state charge transport test.
		double Steady_carrier_density;

		//! Defines the number of events to execute during the equilibration phase of the steady charge transport test.
		int N_equilibration_events;

		//! Specifies whether or not to output the density of occupied states and density of states data after the 
		//! steady state charge transport test.
		bool Enable_state_data_output = true;

		// Exciton Parameters ------------------------------------------------------------------------------------
		// These parameters define the properties of the excitons used by the various simulation tests.

		//! Defines the exciton generation rate on the donor sites is units of cm^-3 s^-3.
		double Exciton_generation_rate_donor; 

		//! Defines the exciton generation rate on the acceptor sites is units of cm^-3 s^-3.
		double Exciton_generation_rate_acceptor;

		//! Defines the singlet exciton lifetime when on the donor sites in units of seconds.
		double Singlet_lifetime_donor; 

		//! Defines the singlet exciton lifetime when on the acceptor sites in units of seconds.
		double Singlet_lifetime_acceptor;

		//! Defines the triplet exciton lifetime when on the donor sites in units of seconds.
		double Triplet_lifetime_donor;

		//! Defines the singlet exciton lifetime when on the acceptor sites in units of seconds.
		double Triplet_lifetime_acceptor;

		//! Defines the singlet exciton hopping rate prefactor when on the donor sites in units of s^-1.
		double R_singlet_hopping_donor;
		
		//! Defines the singlet exciton hopping rate prefactor when on the acceptor sites in units of s^-1.
		double R_singlet_hopping_acceptor;

		//! Defines the singlet exciton localization parameter when on the donor sites in units of nm^-1.
		double Singlet_localization_donor;

		//! Defines the singlet exciton localization parameter when on the acceptor sites in units of nm^-1.
		double Singlet_localization_acceptor;

		//! Defines the triplet exciton hopping rate prefactor when on the donor sites in units of s^-1.
		double R_triplet_hopping_donor;
		
		//! Defines the triplet exciton hopping rate prefactor when on the acceptor sites in units of s^-1.
		double R_triplet_hopping_acceptor;

		//! Defines the triplet exciton localization parameter when on the donor sites in units of nm^-1.
		double Triplet_localization_donor;
		
		//! Defines the triplet exciton localization parameter when on the acceptor sites in units of nm^-1.
		double Triplet_localization_acceptor;

		//! Specifies whether or not to use the FRET model for triplet annihilation events.
		bool Enable_FRET_triplet_annihilation;

		//! Defines the exciton-exciton annihilation rate prefactor when starting from a donor site in units of s^-1.
		double R_exciton_exciton_annihilation_donor;

		//! Defines the exciton-exciton annihilation rate prefactor when starting from an acceptor site in units of s^-1.
		double R_exciton_exciton_annihilation_acceptor;

		//! Defines the exciton-polaron annihilation rate prefactor when starting from a donor site in units of s^-1.
		double R_exciton_polaron_annihilation_donor;

		//! Defines the exciton-polaron annihilation rate prefactor when starting from an acceptor site in units of s^-1.
		double R_exciton_polaron_annihilation_acceptor;

		//! Defines the cutoff radius for FRET-based exciton mechanisms in units of nm.
		int FRET_cutoff;

		//! Defines the binding energy of excitons on donor sites in units of eV.
		double E_exciton_binding_donor;

		//! Defines the binding energy of excitons on acceptor sites in units of eV.
		double E_exciton_binding_acceptor;

		//! Defines the exciton dissociation rate prefactor for excitons starting on a donor site in units of s^-1.
		double R_exciton_dissociation_donor;

		//! Defines the exciton dissociation rate prefactor for excitons starting on an acceptor site in units of s^-1.
		double R_exciton_dissociation_acceptor;

		//! Defines the cutoff radius for exciton dissociation events in units of nm.
		int Exciton_dissociation_cutoff;

		//! Defines the intersystem crossing (singlet to triplet) rate constant of excitons on donor sites in units of s^-1.
		double R_exciton_isc_donor;
		
		//! Defines the intersystem crossing (singlet to triplet) rate constant of excitons on acceptor sites in units of s^-1.
		double R_exciton_isc_acceptor;
		
		//! Defines the reverse intersystem crossing (triplet to singlet) rate constant of excitons on donor sites in units of s^-1.
		double R_exciton_risc_donor;

		//! Defines the reverse intersystem crossing (triplet to singlet) rate constant of excitons on acceptor sites in units of s^-1.
		double R_exciton_risc_acceptor;

		//! Defines the potential energy difference between singlet and triplet exciton states on donor sites in units of eV.
		double E_exciton_ST_donor;
		
		//! Defines the potential energy difference between singlet and triplet exciton states on acceptor sites in units of eV.
		double E_exciton_ST_acceptor;

		// Polaron Parameters ------------------------------------------------------------------------------------------
		// These parameters define the properties of the polarons used by the various simulation tests.

		//! Specifies whether or not holes will be restricted to donor sites and electrons to acceptor sites. 
		bool Enable_phase_restriction;

		//! Defines the polaron hopping rate prefactor of polarons on donor sites in units of s^-1.
		double R_polaron_hopping_donor;
		
		//! Defines the polaron hopping rate prefactor of polarons on acceptor sites in units of s^-1.
		double R_polaron_hopping_acceptor;

		//! Defines the polaron localization parameter when on the donor sites in units of nm^-1.
		double Polaron_localization_donor;

		//! Defines the polaron localization parameter when on the acceptor sites in units of nm^-1.
		double Polaron_localization_acceptor;

		//! Specifies whether or not to use the Miller-Abrahams model for polaron hopping events.
		bool Enable_miller_abrahams;

		//! Specifies whether or not the use the Marcus model for polaron hopping events.
		bool Enable_marcus;

		//! Defines the reorganization energy parameter used by the Marcus model for polaron hopping rates when on the donor sites in units of eV.
		double Reorganization_donor;

		//! Defines the reorganization energy parameter used by the Marcus model for polaron hopping rates when on the acceptor sites in units of eV.
		double Reorganization_acceptor;

		//! Defines the polaron recombination rate prefactor in units of s^-1.
		double R_polaron_recombination;

		//! Defines the cutoff radius of all polaron hopping events in units of nm.
		int Polaron_hopping_cutoff;

		//! Specifies whether or not to use the Gaussian polaron delocalization model.
		bool Enable_gaussian_polaron_delocalization;

		//! Defines the polaron delocalization length used by the Gaussian polaron delocalization model.
		double Polaron_delocalization_length;

		// Additional Lattice Parameters -------------------------------------------------------------------------
		// These parameters define the properties of the sites in the lattice.

		//! Defines the highest occupied molecular orbital energy of the donor sites in units of eV.
		double Homo_donor;

		//! Defines the lowest unoccupied molecular orbital energy of the donor sites in units of eV.
		double Lumo_donor;

		//! Defines the highest occupied molecular orbital energy of the acceptor sites in units of eV.
		double Homo_acceptor;

		//! Defines the lowest unoccupied molecular orbital energy of the acceptor sites in units of eV.
		double Lumo_acceptor;

		//! Specifies whether or not to use the Gaussian density of states model to define the site energies.
		bool Enable_gaussian_dos;

		//! Defines the standard deviation of the donor site energy distribution used by the Gaussian density of states model in units of eV.
		double Energy_stdev_donor;

		//! Defines the standard deviation of the acceptor site energy distribution used by the Gaussian density of states model in units of eV.
		double Energy_stdev_acceptor;

		//! Specifies whether or not to use the exponential density of states model to define the site energies.
		bool Enable_exponential_dos;

		//! Defines the Urbach energy of the donor site energy distribution used by the exponential density of states model in units of eV.
		double Energy_urbach_donor;

		//! Defines the Urbach energy of the acceptor site energy distribution used by the exponential density of states model in units of eV.
		double Energy_urbach_acceptor;

		//! Specifies whether or not the use the correlated Gaussian disorder model to define the site energies.
		bool Enable_correlated_disorder;

		//! Defines the desired correlation length for the correlated Gaussian disorder model in units of nm.
		double Disorder_correlation_length;

		//! Specifies whether or not to use the Gaussian kernel with the correlated Gaussian density of states model.
		bool Enable_gaussian_kernel;

		//! Specifies whether or not to use the power kernel with the correlated Gaussian density of states model.
		bool Enable_power_kernel;

		//! Defines the power kernel exponent to be used with the power kernel correlated Gaussian density of states model.
		int Power_kernel_exponent;

		//! Specifies whether or not to use the interfacial energy shift model to modify the site energies.
		bool Enable_interfacial_energy_shift;

		//! Defines the energy shift parameter for the donor sites used by the interfacial energy shift model in units of eV.
		double Energy_shift_donor;

		//! Defines the energy shift parameter for the acceptor sites used by the interfacial energy shift model in units of eV.
		double Energy_shift_acceptor;

		//! Specifies whether or not to import the site energies from a text file.
		bool Enable_import_energies;

		//! Specifies whether or not to import the site occupancies from a text file.
		bool Enable_import_occupancies;

		//! Specifies whether or not to export the site energies from a text file.
		bool Enable_export_energies;

        //! Specifies whether or not to export the site occupancies from a text file.
		bool Enable_export_occupancies;

        //! Specifies whether only the newest occupancie file, or all of them are kept.
        bool Keep_only_newest_occupancy;

        //! Specifies the interval between two occupancie exports in number of events.
        int Output_interval;

		//! The name of site energy text file to be imported.
		std::string Energies_import_filename;

		//! The name of site occupancies text file to be imported.
		std::string Occupancies_import_filename;

		//! The name of site energy text file to be exported.
		std::string Energies_export_filename;

		//! The name of site occupancies text file to be exported.
		std::string Occupancies_export_filename;

		// Coulomb Calculation Parameters ------------------------------------------------------------------------------

		//! Defines the dielectric constant of the donor sites.
		double Dielectric_donor;

		//! Defines the dielectric constant of the acceptor sites.
		double Dielectric_acceptor;

		//! Defines the cutoff radius for Coulomb interactions in units of nm.
		int Coulomb_cutoff;

		// Functions ---------------------------------------------------------------------------------------------------

		//! \brief Checks the validity of the current parameter values.
		//! \return true if all of the parameter values are valid.
		//! \return false if any of the parameter values are invalid.
		bool checkParameters() const;

		//! \brief Imports the values for all parameters by parsing the input filestream.
		//! \param inputfile is an opened input filestream pointing to a valid parameter text file.
		//! \return true if the parameters are successfully imported.
		//! \return false if there is an error trying to import the parameters.
		bool importParameters(std::ifstream& inputfile);

	private:

	};

}

#endif // EXCIMONTEC_PARAMETERS_H
