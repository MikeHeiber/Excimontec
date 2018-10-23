// Copyright (c) 2017-2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "gtest/gtest.h"
#include "OSC_Sim.h"
#include "Exciton.h"
#include "Utils.h"

using namespace std;
using namespace KMC_Lattice;
using namespace Excimontec;

namespace OSC_SimTests {

	class OSC_SimTest : public ::testing::Test {
	protected:
		Parameters_OPV params_default;
		Parameters_Simulation params_base;
		OSC_Sim sim;
		void SetUp() {
			// General Parameters
			params_default.Enable_FRM = false;
			params_default.Enable_selective_recalc = true;
			params_default.Recalc_cutoff = 1;
			params_default.Enable_full_recalc = false;
			params_default.Enable_logging = false;
			params_default.Enable_periodic_x = true;
			params_default.Enable_periodic_y = true;
			params_default.Enable_periodic_z = true;
			params_default.Length = 100;
			params_default.Width = 100;
			params_default.Height = 100;
			params_default.Unit_size = 1.0;
			params_default.Temperature = 300;
			// Output files
			params_default.Logfile = NULL;
			// Additional General Parameters
			params_default.Internal_potential = 0.0;
			// Morphology Parameters
			params_default.Enable_neat = true;
			params_default.Enable_bilayer = false;
			params_default.Thickness_donor = 50;
			params_default.Thickness_acceptor = 50;
			params_default.Enable_random_blend = false;
			params_default.Acceptor_conc = 0.5;
			params_default.Enable_import_morphology = false;
			params_default.Morphology_filename = "";
			// Test Parameters
			params_default.N_tests = 1000;
			params_default.Enable_exciton_diffusion_test = true;
			params_default.Enable_ToF_test = false;
			params_default.ToF_polaron_type = true;
			params_default.ToF_initial_polarons = 5;
			params_default.Enable_ToF_random_placement = true;
			params_default.Enable_ToF_energy_placement = false;
			params_default.ToF_transient_start = 1e-10;
			params_default.ToF_transient_end = 1e-4;
			params_default.ToF_pnts_per_decade = 20;
			params_default.Enable_IQE_test = false;
			params_default.IQE_time_cutoff = 1e-4;
			params_default.Enable_dynamics_test = false;
			params_default.Enable_dynamics_extraction = false;
			params_default.Dynamics_initial_exciton_conc = 1e16;
			params_default.Dynamics_transient_start = 1e-13;
			params_default.Dynamics_transient_end = 1e-5;
			params_default.Dynamics_pnts_per_decade = 20;
			// Exciton Parameters
			params_default.Exciton_generation_rate_donor = 1e18;
			params_default.Exciton_generation_rate_acceptor = 1e18;
			params_default.Singlet_lifetime_donor = 1e-9;
			params_default.Singlet_lifetime_acceptor = 1e-9;
			params_default.Triplet_lifetime_donor = 1e-8;
			params_default.Triplet_lifetime_acceptor = 1e-8;
			params_default.R_singlet_hopping_donor = 1e11;
			params_default.R_singlet_hopping_acceptor = 1e11;
			params_default.Singlet_localization_donor = 1.0;
			params_default.Singlet_localization_acceptor = 1.0;
			params_default.R_triplet_hopping_donor = 1e11;
			params_default.R_triplet_hopping_acceptor = 1e11;
			params_default.Triplet_localization_donor = 1.0;
			params_default.Triplet_localization_acceptor = 1.0;
			params_default.Enable_FRET_triplet_annihilation = false;
			params_default.R_exciton_exciton_annihilation_donor = 1e1;
			params_default.R_exciton_exciton_annihilation_acceptor = 1e1;
			params_default.R_exciton_polaron_annihilation_donor = 1e1;
			params_default.R_exciton_polaron_annihilation_acceptor = 1e1;
			params_default.FRET_cutoff = 1;
			params_default.E_exciton_binding_donor = 0.5;
			params_default.E_exciton_binding_acceptor = 0.5;
			params_default.R_exciton_dissociation_donor = 1e14;
			params_default.R_exciton_dissociation_acceptor = 1e14;
			params_default.Exciton_dissociation_cutoff = 1;
			params_default.R_exciton_isc_donor = 1e-12;
			params_default.R_exciton_isc_acceptor = 1e-12;
			params_default.R_exciton_risc_donor = 1e-12;
			params_default.R_exciton_risc_acceptor = 1e-12;
			params_default.E_exciton_ST_donor = 0.7;
			params_default.E_exciton_ST_acceptor = 0.7;
			// Polaron Parameters
			params_default.Enable_phase_restriction = true;
			params_default.R_polaron_hopping_donor = 1e11;
			params_default.R_polaron_hopping_acceptor = 1e11;
			params_default.Polaron_localization_donor = 2.0;
			params_default.Polaron_localization_acceptor = 2.0;
			params_default.Enable_miller_abrahams = true;
			params_default.Enable_marcus = false;
			params_default.Reorganization_donor = 0.2;
			params_default.Reorganization_acceptor = 0.2;
			params_default.R_polaron_recombination = 1e10;
			params_default.Polaron_hopping_cutoff = 1;
			params_default.Enable_gaussian_polaron_delocalization = false;
			params_default.Polaron_delocalization_length = 1.0;
			// Additional Lattice Parameters
			params_default.Homo_donor = 5.0;
			params_default.Lumo_donor = 3.0;
			params_default.Homo_acceptor = 6.0;
			params_default.Lumo_acceptor = 4.0;
			params_default.Enable_gaussian_dos = false;
			params_default.Energy_stdev_donor = 0.05;
			params_default.Energy_stdev_acceptor = 0.05;
			params_default.Enable_exponential_dos = false;
			params_default.Energy_urbach_donor = 0.03;
			params_default.Energy_urbach_acceptor = 0.03;
			params_default.Enable_correlated_disorder = false;
			params_default.Disorder_correlation_length = 1.0;
			params_default.Enable_gaussian_kernel = false;
			params_default.Enable_power_kernel = false;
			params_default.Power_kernel_exponent = -2;
			// Coulomb Calculation Parameters
			params_default.Dielectric_donor = 3.5;
			params_default.Dielectric_acceptor = 3.5;
			params_default.Coulomb_cutoff = 50;
		}
	};

	TEST_F(OSC_SimTest, ParameterCheckTests) {
		sim = OSC_Sim();
		// Check that default parameters are valid
		EXPECT_TRUE(sim.init(params_default, 0));
		// Check various invalid parameter sets
		Parameters_OPV params;
		// Check that a device architecture is defined
		params.Enable_neat = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for multiple enabled morphologies
		params = params_default;
		params.Enable_neat = true;
		params.Enable_random_blend = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for valid bilayer thicknesses
		params = params_default;
		params.Enable_neat = false;
		params.Enable_bilayer = true;
		params.Thickness_donor = 20;
		params.Thickness_acceptor = 20;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for invalid lattice dimensions
		params = params_default;
		params.Height = 0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for invalid unit size
		params = params_default;
		params.Unit_size = 0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for invalid temp
		params = params_default;
		params.Temperature = 0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for valid selection of KMC algorithm
		params = params_default;
		// Check for multiple KMC algorithms enabled
		params.Enable_FRM = true;
		params.Enable_selective_recalc = true;
		params.Enable_full_recalc = false;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_FRM = true;
		params.Enable_selective_recalc = false;
		params.Enable_full_recalc = true;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_FRM = false;
		params.Enable_selective_recalc = true;
		params.Enable_full_recalc = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for no KMC algorithm specified
		params.Enable_FRM = false;
		params.Enable_selective_recalc = false;
		params.Enable_full_recalc = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check recalculation radius when using the selective recalc method
		params.Enable_selective_recalc = true;
		params.Recalc_cutoff = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params.Recalc_cutoff = 1;
		params.FRET_cutoff = 2;
		params.Exciton_dissociation_cutoff = 1;
		params.Polaron_hopping_cutoff = 1;
		EXPECT_FALSE(sim.init(params, 0));
		params.Recalc_cutoff = 1;
		params.FRET_cutoff = 1;
		params.Exciton_dissociation_cutoff = 2;
		params.Polaron_hopping_cutoff = 1;
		EXPECT_FALSE(sim.init(params, 0));
		params.Recalc_cutoff = 1;
		params.FRET_cutoff = 1;
		params.Exciton_dissociation_cutoff = 1;
		params.Polaron_hopping_cutoff = 2;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for multiple enabled tests
		params = params_default;
		params.Enable_dynamics_test = true;
		params.Enable_exciton_diffusion_test = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for no tests enabled
		params.Enable_dynamics_test = false;
		params.Enable_exciton_diffusion_test = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for invalid test duration
		params.Enable_exciton_diffusion_test = true;
		params.N_tests = 0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for valid ToF test options
		// Check valid parameters
		params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_ToF_test = true;
		params.Enable_periodic_z = false;
		EXPECT_TRUE(sim.init(params, 0));
		// Check for invalid z boundary conditions
		params.Enable_periodic_z = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for invalid morphology
		params.Enable_periodic_z = false;
		params.Enable_bilayer = true;
		params.Enable_neat = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for valid placement option
		params.Enable_neat = true;
		params.Enable_bilayer = false;
		params.Enable_ToF_random_placement = false;
		params.Enable_ToF_energy_placement = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for invalid polaron type
		params.Enable_ToF_random_placement = true;
		params.ToF_polaron_type = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for valid IQE test options
		// Check valid parameters
		params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_IQE_test = true;
		params.Enable_periodic_z = false;
		params.Enable_neat = false;
		params.Enable_bilayer = true;
		EXPECT_TRUE(sim.init(params, 0));
		// Check invalid boundary conditions
		params.Enable_periodic_z = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check invalid device architecture
		params.Enable_periodic_z = false;
		params.Enable_neat = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for valid Dynamics test options
		params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_dynamics_test = true;
		// Check for dynamics extraction test and internal potential
		params.Enable_dynamics_extraction = false;
		params.Enable_periodic_z = false;
		params.Internal_potential = -1.0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for dynamics extraction test and z boundary conditions
		params.Enable_dynamics_extraction = false;
		params.Enable_periodic_z = true;
		params.Internal_potential = -1.0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check exciton parameter values
		params = params_default;
		params.Exciton_generation_rate_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Singlet_lifetime_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Triplet_lifetime_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_singlet_hopping_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Singlet_localization_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_triplet_hopping_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Triplet_localization_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_exciton_exciton_annihilation_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_exciton_polaron_annihilation_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.FRET_cutoff = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.E_exciton_binding_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_exciton_dissociation_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Exciton_dissociation_cutoff = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_exciton_isc_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_exciton_risc_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.E_exciton_ST_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		// Check polaron parameter values
		params = params_default;
		params.R_polaron_hopping_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Polaron_localization_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_miller_abrahams = true;
		params.Enable_marcus = true;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_miller_abrahams = false;
		params.Enable_marcus = false;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Reorganization_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.R_polaron_recombination = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Polaron_hopping_cutoff = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Polaron_delocalization_length = 0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check lattice site parameters
		params = params_default;
		params.Homo_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Lumo_acceptor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Enable_exponential_dos = true;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Energy_stdev_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_exponential_dos = true;
		params.Energy_urbach_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_gaussian_dos = false;
		params.Enable_correlated_disorder = true;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Enable_correlated_disorder = true;
		params.Disorder_correlation_length = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Enable_correlated_disorder = true;
		params.Enable_gaussian_kernel = true;
		params.Disorder_correlation_length = 2.5;
		EXPECT_FALSE(sim.init(params, 0));
		// Coulomb calculation options
		params = params_default;
		params.Dielectric_donor = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params = params_default;
		params.Coulomb_cutoff = 0;
		EXPECT_FALSE(sim.init(params, 0));
	}

	TEST_F(OSC_SimTest, SetupTests) {
		sim = OSC_Sim();
		Parameters_OPV params = params_default;
		EXPECT_TRUE(sim.init(params, 0));
		EXPECT_EQ(params.Length*params.Width*params.Height*1e-21, sim.getVolume());
		// Test morphology import
		sim = OSC_Sim();
		params.Enable_neat = false;
		params.Enable_import_morphology = true;
		params.Length = 200;
		params.Width = 200;
		params.Height = 200;
		// Test incorrect filename
		params.Morphology_filename = "./test/test_morphology123.txt";
		EXPECT_FALSE(sim.init(params, 0));
		params.Morphology_filename = "./test/test_morphology.txt";
		// Test incorrect dimensions
		sim = OSC_Sim();
		params.Height = 100;
		EXPECT_FALSE(sim.init(params, 0));
		params.Height = 200;
		// Test correct import
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
	}

	TEST_F(OSC_SimTest, ExcitonDynamicsTests) {
		// Singlet exciton lifetime test
		sim = OSC_Sim();
		Parameters_OPV params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_dynamics_test = true;
		params.N_tests = 3000;
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		auto time_data = sim.getDynamicsTransientTimes();
		auto singlet_data = sim.getDynamicsTransientSinglets();
		EXPECT_TRUE(time_data.size() == singlet_data.size());
		vector<pair<double, double>> transient_data(time_data.size());
		for (int i = 0; i < (int)time_data.size(); i++) {
			transient_data[i] = make_pair(time_data[i], (double)singlet_data[i] / (sim.getVolume()*sim.getN_transient_cycles()));
		}
		EXPECT_NEAR(params.Dynamics_initial_exciton_conc, transient_data[0].second, 1e-3*params.Dynamics_initial_exciton_conc);
		EXPECT_NEAR(1 / exp(1), interpolateData(transient_data, params.Singlet_lifetime_donor) / params.Dynamics_initial_exciton_conc, 3e-2);
		// Triplet exciton lifetime test
		sim = OSC_Sim();
		params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_dynamics_test = true;
		params.R_exciton_isc_donor = 1e16;
		params.N_tests = 3000;
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		time_data = sim.getDynamicsTransientTimes();
		auto triplet_data = sim.getDynamicsTransientTriplets();
		EXPECT_TRUE(time_data.size() == triplet_data.size());
		transient_data.resize(time_data.size());
		for (int i = 0; i < (int)time_data.size(); i++) {
			transient_data[i] = make_pair(time_data[i], (double)triplet_data[i] / (sim.getVolume()*sim.getN_transient_cycles()));
		}
		EXPECT_NEAR(1 / exp(1), interpolateData(transient_data, params.Triplet_lifetime_donor) / params.Dynamics_initial_exciton_conc, 3e-2);
	}

	TEST_F(OSC_SimTest, ExcitonDiffusionTests) {
		// Singlet exciton diffusion test
		sim = OSC_Sim();
		Parameters_OPV params = params_default;
		params.N_tests = 5000;
		EXPECT_TRUE(sim.init(params, 0));
		bool success;
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		// Check exciton creation statistics
		EXPECT_EQ(params.N_tests, sim.getN_excitons_created());
		EXPECT_EQ(params.N_tests, sim.getN_excitons_created((char)1));
		EXPECT_EQ(0, sim.getN_excitons_created((char)2));
		// Check exciton diffusion results
		auto lifetime_data = sim.getExcitonLifetimeData();
		EXPECT_NEAR(params.Singlet_lifetime_donor, vector_avg(lifetime_data), 5e-2*params.Singlet_lifetime_donor);
		EXPECT_DOUBLE_EQ(1.0, vector_avg(sim.getExcitonHopLengthData()));
		auto displacement_data = sim.getExcitonDiffusionData();
		auto ratio_data = displacement_data;
		transform(displacement_data.begin(), displacement_data.end(), lifetime_data.begin(), ratio_data.begin(), [params](double& displacement_element, double& lifetime_element) {
			return displacement_element / sqrt(6 * params.R_singlet_hopping_donor*lifetime_element);
		});
		double dim = 3.0;
		double expected_ratio = sqrt(2.0 / dim)*(tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0));
		EXPECT_NEAR(expected_ratio, vector_avg(ratio_data), 2e-2*expected_ratio);
		// Triplet exciton diffusion test
		sim = OSC_Sim();
		params = params_default;
		params.R_exciton_isc_donor = 1e16;
		params.N_tests = 5000;
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		lifetime_data = sim.getExcitonLifetimeData();
		EXPECT_NEAR(params.Triplet_lifetime_donor, vector_avg(lifetime_data), 5e-2*params.Triplet_lifetime_donor);
		EXPECT_DOUBLE_EQ(1.0, vector_avg(sim.getExcitonHopLengthData()));
		displacement_data = sim.getExcitonDiffusionData();
		ratio_data = displacement_data;
		transform(displacement_data.begin(), displacement_data.end(), lifetime_data.begin(), ratio_data.begin(), [params](double& displacement_element, double& lifetime_element) {
			return displacement_element / sqrt(6 * params.R_triplet_hopping_donor * exp(-2.0 * params.Triplet_localization_donor) * lifetime_element);
		});
		dim = 3.0;
		expected_ratio = sqrt(2.0 / dim)*(tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0));
		EXPECT_NEAR(expected_ratio, vector_avg(ratio_data), 2e-2*expected_ratio);
	}

	TEST_F(OSC_SimTest, IQETests) {
		// Setup starting parameters
		Parameters_OPV params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_IQE_test = true;
		params.Enable_neat = false;
		params.Enable_periodic_z = false;
		params.Enable_bilayer = true;
		params.Height = 40;
		params.Thickness_donor = 20;
		params.Thickness_acceptor = 20;
		params.FRET_cutoff = 4;
		params.Polaron_hopping_cutoff = 3;
		params.Recalc_cutoff = 4;
		params.Internal_potential = -1.0;
		params.N_tests = 500;
		bool success;
		// Baseline bilayer IQE test
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		double IQE1 = 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2.0 * (double)sim.getN_excitons_created());
		// Check for field activated charge separation
		params.Internal_potential = -2.0;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		double IQE2 = 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2.0 * (double)sim.getN_excitons_created());
		EXPECT_GT(IQE2, IQE1);
		// Check for thermally activated charge separation
		params.Internal_potential = -1.0;
		params.Temperature = 350;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		double IQE3 = 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2.0 * (double)sim.getN_excitons_created());
		EXPECT_GT(IQE3, IQE1);
		// Check for delocalization enhanced charge separation
		params.Temperature = 300;
		params.Enable_gaussian_polaron_delocalization = true;
		params.Polaron_delocalization_length = 2.0;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		double IQE4 = 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2.0 * (double)sim.getN_excitons_created());
		EXPECT_GT(IQE4, IQE1);
		// Check for recombination rate dependent charge separation
		params.Enable_gaussian_polaron_delocalization = false;
		params.R_polaron_recombination = 1e9;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		double IQE5 = 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2.0 * (double)sim.getN_excitons_created());
		EXPECT_GT(IQE5, IQE1);
	}

	TEST_F(OSC_SimTest, ToFTests) {
		// Hole ToF test
		sim = OSC_Sim();
		Parameters_OPV params = params_default;
		params.Enable_periodic_z = false;
		params.Height = 200;
		params.Internal_potential = -4.0;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_ToF_test = true;
		params.N_tests = 1000;
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		auto transit_time_data = sim.getTransitTimeData();
		// Check that the transit time probability histogram sums to 1
		auto hist = sim.calculateTransitTimeHist(transit_time_data, (int)transit_time_data.size());
		auto cum_hist = calculateCumulativeHist(hist);
		EXPECT_NEAR(1.0, cum_hist.back().second, 1e-3);
		// Check the mobility compared to analytical expectation
		auto mobility_data = sim.calculateMobilityData(transit_time_data);
		double dim = 3.0;
		double expected_mobility = (params.R_polaron_hopping_donor*exp(-2.0*params.Polaron_localization_donor)*1e-14) * (2.0 / 3.0) * (tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0)) * (1 / (K_b*params.Temperature));
		EXPECT_NEAR(expected_mobility, vector_avg(mobility_data), 1e-1*expected_mobility);
		// Electron ToF test on neat should not be allowed
		sim = OSC_Sim();
		params = params_default;
		params.Enable_periodic_z = false;
		params.Height = 200;
		params.Internal_potential = -4.0;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_ToF_test = true;
		params.ToF_polaron_type = false;
		params.N_tests = 1000;
		EXPECT_FALSE(sim.init(params, 0));
		// Electron ToF test on random blend should work
		sim = OSC_Sim();
		params.Enable_neat = false;
		params.Enable_random_blend = true;
		params.Acceptor_conc = 0.99; // Dilute blend should behave like a neat electron transport material might
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		transit_time_data = sim.getTransitTimeData();
		mobility_data = sim.calculateMobilityData(transit_time_data);
		dim = 3.0;
		expected_mobility = (params.R_polaron_hopping_donor*exp(-2.0*params.Polaron_localization_donor)*1e-14) * (2.0 / 3.0) * (tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0)) * (1 / (K_b*params.Temperature));
		EXPECT_NEAR(expected_mobility, vector_avg(mobility_data), 1e-1*expected_mobility);
	}

	TEST_F(OSC_SimTest, CorrelatedDisorderGaussianKernelTests) {
		sim = OSC_Sim();
		Parameters_OPV params = params_default;
		params.Enable_gaussian_dos = true;
		params.Energy_stdev_donor = 0.05;
		params.Energy_stdev_acceptor = 0.05;
		params.Enable_correlated_disorder = true;
		params.Enable_gaussian_kernel = true;
		params.Enable_power_kernel = false;
		// correlation length = 1.0, Unit_size = 0.8
		params.Length = 40;
		params.Width = 40;
		params.Height = 40;
		params.Disorder_correlation_length = 1.1;
		sim.init(params, 0);
		auto energies = sim.getSiteEnergies(1);
		EXPECT_NEAR(0.0, vector_avg(energies), 5e-3);
		EXPECT_NEAR(0.05, vector_stdev(energies), 1e-3);
		auto correlation_data = sim.getDOSCorrelationData();
		EXPECT_NEAR(1 / exp(1), interpolateData(correlation_data, 1.1), 0.05);
		// correlation length = 1.3, Unit_size = 1.2
		sim = OSC_Sim();
		params.Unit_size = 1.2;
		params.Length = 40;
		params.Width = 40;
		params.Height = 40;
		params.Disorder_correlation_length = 1.3;
		sim.init(params, 0);
		energies = sim.getSiteEnergies(1);
		EXPECT_NEAR(0.0, vector_avg(energies), 5e-3);
		EXPECT_NEAR(0.05, vector_stdev(energies), 1e-3);
		correlation_data = sim.getDOSCorrelationData();
		EXPECT_NEAR(1 / exp(1), interpolateData(correlation_data, 1.3), 0.05);
	}

}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	// Redirect cout to NULL to suppress command line output during the tests
	//cout.rdbuf(NULL);
	return RUN_ALL_TESTS();
}
