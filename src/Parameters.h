// Copyright (c) 2017-2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Lattice.h"
#include "Simulation.h"
#include "Utils.h"



namespace Excimontec {

	//! \brief This class extends the Object class to create an exciton object to represent a singlet or triplet exciton in an organic semiconductor.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2017-2018
	class Parameters : public KMC_Lattice::Parameters_Simulation {
	public:
		// main parameters
		bool Enable_import_morphology_single;
		bool Enable_import_morphology_set;
		std::string Morphology_set_format;
		int N_test_morphologies;
		int N_morphology_set_size;
		bool Enable_extraction_map_output;
		// Additional General Parameters
		double Internal_potential;
		// Morphology Parameters
		bool Enable_neat; // Neat takes on donor properties
		bool Enable_bilayer;
		int Thickness_donor; // sites
		int Thickness_acceptor; // sites
		bool Enable_random_blend;
		double Acceptor_conc;
		bool Enable_import_morphology;
		std::string Morphology_filename;
		// Test Parameters
		int N_tests;
		bool Enable_exciton_diffusion_test;
		bool Enable_ToF_test;
		bool ToF_polaron_type;
		int ToF_initial_polarons;
		bool Enable_ToF_random_placement;
		bool Enable_ToF_energy_placement;
		double ToF_placement_energy;
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
		double Singlet_localization_donor;
		double Singlet_localization_acceptor;
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
		bool Enable_interfacial_energy_shift;
		double Energy_shift_donor;
		double Energy_shift_acceptor;
		bool Enable_import_energies;
		std::string Energies_import_filename;
		// Coulomb Calculation Parameters
		double Dielectric_donor;
		double Dielectric_acceptor;
		int Coulomb_cutoff; // nm

		// Functions
		bool checkParameters() const;
		bool importParameters(std::ifstream& inputfile);

	private:
	
	};

}

#endif // PARAMETERS_H
