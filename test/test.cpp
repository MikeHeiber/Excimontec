// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "gtest/gtest.h"
#include "OSC_Sim.h"
#include "Exciton.h"
#include "Utils.h"

using namespace std;
using namespace Utils;

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
			params_default.Recalc_cutoff = 3;
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
			params_default.Morphology_file = NULL;
			// Test Parameters
			params_default.N_tests = 10000;
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
			params_default.Singlet_lifetime_donor = 500e-12;
			params_default.Singlet_lifetime_acceptor = 500e-12;
			params_default.Triplet_lifetime_donor = 1e-6;
			params_default.Triplet_lifetime_acceptor = 1e-6;
			params_default.R_singlet_hopping_donor = 1e11;
			params_default.R_singlet_hopping_acceptor = 1e11;
			params_default.Singlet_localization_donor = 1.0;
			params_default.Singlet_localization_acceptor = 1.0;
			params_default.R_singlet_hopping_acceptor = 1e12;
			params_default.R_triplet_hopping_donor = 1e12;
			params_default.R_triplet_hopping_acceptor = 1e12;
			params_default.Triplet_localization_donor = 2.0;
			params_default.Triplet_localization_acceptor = 2.0;
			params_default.Enable_FRET_triplet_annihilation = false;
			params_default.R_exciton_exciton_annihilation_donor = 1e12;
			params_default.R_exciton_exciton_annihilation_acceptor = 1e12;
			params_default.R_exciton_polaron_annihilation_donor = 1e12;
			params_default.R_exciton_polaron_annihilation_acceptor = 1e12;
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
			params_default.R_polaron_recombination = 1e12;
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
		Parameters_OPV params = params_default;
		// Check that default parameters are valid
		EXPECT_TRUE(sim.init(params, 0));
		// Check for multiple enabled morphologies
		params.Enable_neat = true;
		params.Enable_random_blend = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for ToF test boundary conditions, films architectures, and test conditions
		sim = OSC_Sim();
		params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_ToF_test = true;
		params.Enable_periodic_z = false;
		EXPECT_TRUE(sim.init(params, 0));
		params.Enable_periodic_z = true;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_periodic_z = false;
		params.Enable_bilayer = true;
		params.Enable_neat = false;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_neat = true;
		params.Enable_bilayer = false;
		params.Enable_ToF_random_placement = false;
		params.Enable_ToF_energy_placement = false;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_ToF_random_placement = true;
		params.ToF_polaron_type = false;
		EXPECT_FALSE(sim.init(params, 0));
		// Check for multiple enabled tests
		params = params_default;
		params.Enable_dynamics_test = true;
		params.Enable_exciton_diffusion_test = true;
		EXPECT_FALSE(sim.init(params, 0));
	}

	TEST_F(OSC_SimTest, ExcitonDiffusionTest) {
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params_default, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		auto lifetime_data = sim.getExcitonLifetimeData();
		double lifetime_avg = vector_avg(sim.getExcitonLifetimeData());
		EXPECT_NEAR(params_default.Singlet_lifetime_donor, lifetime_avg, 2e-2*params_default.Singlet_lifetime_donor);
		double hop_length_avg = vector_avg(sim.getExcitonHopLengthData());
		EXPECT_DOUBLE_EQ(1.0, hop_length_avg);
		auto displacement_data = sim.getExcitonDiffusionData();
		auto ratio_data(displacement_data);
		transform(displacement_data.begin(), displacement_data.end(), lifetime_data.begin(), ratio_data.begin(), [this](double& displacement_element, double& lifetime_element) {
			return displacement_element / sqrt(6 * params_default.R_singlet_hopping_donor*lifetime_element);
		});
		double dim = 3.0;
		double expected_ratio = sqrt(2.0 / dim)*(tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0));
		EXPECT_NEAR(expected_ratio, vector_avg(ratio_data), 2e-2*expected_ratio);
	}

	TEST_F(OSC_SimTest, ToFTest) {
		// Hole ToF test
		sim = OSC_Sim();
		params_default.Enable_periodic_z = false;
		params_default.Height = 200;
		params_default.Internal_potential = -4.0;
		params_default.Enable_exciton_diffusion_test = false;
		params_default.Enable_ToF_test = true;
		params_default.N_tests = 1000;
		EXPECT_TRUE(sim.init(params_default, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		auto transit_time_data = sim.getTransitTimeData();
		auto mobility_data = sim.calculateMobilityData(transit_time_data);
		double dim = 3.0;
		double expected_mobility = (params_default.R_polaron_hopping_donor*exp(-2.0*params_default.Polaron_localization_donor)*1e-14) * (2.0/3.0) * (tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0)) * ( 1 / (K_b*params_default.Temperature));
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
		EXPECT_NEAR(1/exp(1), interpolateData(correlation_data, 1.1), 0.05);
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
	return RUN_ALL_TESTS();
}
