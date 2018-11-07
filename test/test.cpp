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
			params_default.IQE_time_cutoff = 1e-3;
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
			params_default.Enable_interfacial_energy_shift = false;
			params_default.Energy_shift_donor = 0.0;
			params_default.Energy_shift_acceptor = 0.0;
			params_default.Enable_import_energies = false;
			params_default.Energies_import_filename = "energies.txt";
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
		auto params = params_default;
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
		params.Enable_bilayer = false;
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
		params.Enable_dynamics_extraction = true;
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
		// Check enabled correlated disorder without specifying a DOS model
		params = params_default;
		params.Enable_gaussian_dos = false;
		params.Enable_correlated_disorder = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check enabled correlated disorder with exponential DOS model
		params = params_default;
		params.Enable_gaussian_dos = false;
		params.Enable_exponential_dos = true;
		params.Enable_correlated_disorder = true;
		EXPECT_FALSE(sim.init(params, 0));
		// Check correlated disorder correlation length minimum
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Enable_correlated_disorder = true;
		params.Disorder_correlation_length = 0;
		EXPECT_FALSE(sim.init(params, 0));
		// Check correlated disorder correlation length maximum
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Enable_correlated_disorder = true;
		params.Enable_gaussian_kernel = true;
		params.Disorder_correlation_length = 2.5;
		EXPECT_FALSE(sim.init(params, 0));
		// Check correlated disorder with different disorders
		params = params_default;
		params.Enable_neat = false;
		params.Enable_bilayer = true;
		params.Enable_gaussian_dos = true;
		params.Enable_correlated_disorder = true;
		params.Enable_gaussian_kernel = true;
		params.Disorder_correlation_length = 1.0;
		params.Energy_stdev_donor = 0.05;
		params.Energy_stdev_acceptor = 0.025;
		EXPECT_FALSE(sim.init(params, 0));
		// Check correlated disorder power kernel exponent value
		params = params_default;
		params.Enable_gaussian_dos = true;
		params.Enable_correlated_disorder = true;
		params.Enable_power_kernel = true;
		params.Power_kernel_exponent = 0;
		EXPECT_FALSE(sim.init(params, 0));
		params.Power_kernel_exponent = -3;
		EXPECT_FALSE(sim.init(params, 0));
		// Check invalid interfacial energy shift params
		params = params_default;
		params.Enable_interfacial_energy_shift = true;
		params.Energy_shift_donor = -1;
		EXPECT_FALSE(sim.init(params, 0));
		// Check invalid energy import conditions
		// Check missing filename
		params = params_default;
		params.Enable_import_energies = true;
		params.Energies_import_filename = "";
		EXPECT_FALSE(sim.init(params, 0));
		// Check conflicting options
		params.Energies_import_filename = "energies.txt";
		params.Enable_gaussian_dos = true;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_gaussian_dos = false;
		params.Enable_interfacial_energy_shift = true;
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
		// Check default parameters
		sim = OSC_Sim();
		auto params = params_default;
		EXPECT_TRUE(sim.init(params, 0));
		EXPECT_EQ(params.Length*params.Width*params.Height*1e-21, sim.getVolume());
	}

	TEST_F(OSC_SimTest, ObjectCreationTests) {
		sim = OSC_Sim();
		auto params = params_default;
		params.Enable_neat = false;
		params.Enable_bilayer = true;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_IQE_test = true;
		params.Enable_periodic_z = false;
		EXPECT_TRUE(sim.init(params, 0));
		// Create exciton
		sim.createExciton(Coords(10, 10, 10), true);
		// redirect cout to file
		ofstream outfile("./test/status.txt");
		auto old_buf = cout.rdbuf(outfile.rdbuf());
		sim.outputStatus();
		cout.rdbuf(old_buf);
		outfile.close();
		// Read status file
		ifstream infile("./test/status.txt");
		string line;
		vector<string> lines;
		while (getline(infile, line)) {
			lines.push_back(line);
		}
		infile.close();
		EXPECT_EQ(6, (int)lines.size());
		EXPECT_STREQ(lines[3].c_str(), "0: Exciton 1 is at 10,10,10.");
		// Try executing an event
		sim.calculateAllEvents();
		sim.executeNextEvent();
		// redirect cout to file
		ofstream outfile2("./test/status.txt");
		old_buf = cout.rdbuf(outfile2.rdbuf());
		sim.outputStatus();
		cout.rdbuf(old_buf);
		outfile2.close();
		// Read status file
		ifstream infile2("./test/status.txt");
		lines.clear();
		while (getline(infile2, line)) {
			lines.push_back(line);
		}
		infile2.close();
		EXPECT_EQ(6, (int)lines.size());
		EXPECT_STREQ(lines[1].c_str(), "0: 1 excitons have been created and 1 events have been executed.");
		// Create electron and create hole
		sim.createElectron(Coords(10, 10, 49));
		sim.createHole(Coords(10, 10, 50));
		// redirect cout to file
		ofstream outfile3("./test/status.txt");
		old_buf = cout.rdbuf(outfile3.rdbuf());
		sim.outputStatus();
		cout.rdbuf(old_buf);
		outfile3.close();
		// Read status file
		ifstream infile3("./test/status.txt");
		lines.clear();
		while (getline(infile3, line)) {
			lines.push_back(line);
		}
		infile3.close();
		EXPECT_EQ(8, (int)lines.size());
		EXPECT_EQ(lines[5], "0: Electron 1 is at 10,10,49.");
		EXPECT_EQ(lines[7], "0: Hole 1 is at 10,10,50.");
	}

	TEST_F(OSC_SimTest, LoggingTests) {
		sim = OSC_Sim();
		auto params = params_default;
		// Enable logging
		params.Enable_logging = true;
		ofstream logfile("./test/log.txt");
		params.Logfile = &logfile;
		// Run a simple random blend IQE test to test a variety of event types
		params.Enable_neat = false;
		params.Enable_random_blend = true;
		params.Acceptor_conc = 0.3;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_IQE_test = true;
		params.Enable_periodic_z = false;
		params.Internal_potential = -1.0;
		params.R_exciton_exciton_annihilation_donor = 1e14;
		params.R_exciton_exciton_annihilation_acceptor = 1e14;
		params.R_exciton_polaron_annihilation_donor = 1e14;
		params.R_exciton_polaron_annihilation_acceptor = 1e14;
		params.R_exciton_isc_donor = 1e9;
		params.R_exciton_isc_acceptor = 1e9;
		params.R_exciton_risc_donor = 1e9;
		params.R_exciton_risc_acceptor = 1e9;
		params.E_exciton_ST_donor = 0.1;
		params.E_exciton_ST_acceptor = 0.1;
		EXPECT_TRUE(sim.init(params, 0));
		// Execute 1000 events and save events executed
		vector<string> event_types;
		for (int i = 0; i < 1000; i++) {
			sim.executeNextEvent();
			event_types.push_back(sim.getPreviousEventType());
		}
		logfile.close();
		// Extract the events executed lines from the log
		ifstream infile("./test/log.txt");
		string line;
		vector<string> event_lines;
		int n = 0;
		while (getline(infile, line)) {
			if (line.find("Event " + to_string(n)) != string::npos) {
				event_lines.push_back(line);
				n++;
			}
		}
		infile.close();
		// Check that the extracted event lines from the log match the vector of event types
		EXPECT_EQ(event_lines.size(), event_types.size());
		if (event_lines.size() == event_types.size()) {
			for (int i = 0; i < (int)event_lines.size(); i++) {
				EXPECT_TRUE(event_lines[i].find(event_types[i]) != string::npos);
			}
		}
	}

	TEST_F(OSC_SimTest, EnergiesImportTests) {
		// Create sample energies file
		sim = OSC_Sim();
		auto params = params_default;
		params.Length = 30;
		params.Width = 30;
		params.Height = 30;
		params.Enable_gaussian_dos = true;
		EXPECT_TRUE(sim.init(params, 0));
		double energies_stdev1 = vector_stdev(sim.getSiteEnergies(1));
		sim.exportEnergies("./test/energies.txt");
		// Test valid import operation
		sim = OSC_Sim();
		params.Enable_gaussian_dos = false;
		params.Enable_import_energies = true;
		params.Energies_import_filename = "./test/energies.txt";
		EXPECT_TRUE(sim.init(params, 0));
		// Check imported energies
		auto site_energies = sim.getSiteEnergies(1);
		EXPECT_NEAR(0, vector_avg(site_energies), 5e-3);
		EXPECT_NEAR(energies_stdev1, vector_stdev(site_energies), 1e-4);
		// Test missing energies file
		params.Energies_import_filename = "energies.txt";
		EXPECT_FALSE(sim.init(params, 0));
		// Test energies file with missing data
		params.Energies_import_filename = "./test/energies_missing_data.txt";
		EXPECT_FALSE(sim.init(params, 0));
		// Test energy file with no dimensions
		params.Energies_import_filename = "./test/energies_missing_dims.txt";
		EXPECT_FALSE(sim.init(params, 0));
		// Test energy file with improper dimensions
		params.Length = 50;
		params.Width = 50;
		params.Height = 50;
		params.Energies_import_filename = "./test/energies.txt";
		EXPECT_FALSE(sim.init(params, 0));
	}

	TEST_F(OSC_SimTest, MorphologyImportTests) {
		auto params = params_default;
		// Test morphology import
		sim = OSC_Sim();
		params.Enable_neat = false;
		params.Enable_import_morphology = true;
		params.Length = 50;
		params.Width = 50;
		params.Height = 50;
		// Test incorrect filename
		params.Morphology_filename = "./test/test_morphology123.txt";
		EXPECT_FALSE(sim.init(params, 0));
		// Test incorrect dimensions
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_v4-0_compressed.txt";
		params.Height = 100;
		EXPECT_FALSE(sim.init(params, 0));
		params.Height = 50;
		// Test correct import of v3.2 compressed
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_v3-2_compressed.txt";
		EXPECT_TRUE(sim.init(params, 0));
		// Test correct import of v3.2 uncompressed
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_v3-2_uncompressed.txt";
		EXPECT_TRUE(sim.init(params, 0));
		// Test correct import of v4.0 compressed
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_v4-0_compressed.txt";
		EXPECT_TRUE(sim.init(params, 0));
		// Test correct import of v4.0 uncompressed
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_v4-0_uncompressed.txt";
		EXPECT_TRUE(sim.init(params, 0));
		// Test behavior when there is missing data
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_v4-0_missing_data.txt";
		EXPECT_FALSE(sim.init(params, 0));
		// Test behavior when there is no header line
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_no_version.txt";
		EXPECT_FALSE(sim.init(params, 0));
		// Test behavior when the version is too old
		sim = OSC_Sim();
		params.Morphology_filename = "./test/morphology_old_version.txt";
		EXPECT_FALSE(sim.init(params, 0));
	}

	TEST_F(OSC_SimTest, ExcitonDynamicsTests) {
		// Singlet exciton lifetime test
		sim = OSC_Sim();
		Parameters_OPV params = params_default;
		params.Enable_exciton_diffusion_test = false;
		params.Enable_dynamics_test = true;
		params.Dynamics_transient_end = 1e-6;
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
		params.Length = 50;
		params.Width = 50;
		params.Height = 40;
		params.Coulomb_cutoff = 25;
		params.Thickness_donor = 20;
		params.Thickness_acceptor = 20;
		params.Recalc_cutoff = 2;
		params.FRET_cutoff = 2;
		params.Exciton_dissociation_cutoff = 2;
		params.Polaron_hopping_cutoff = 2;
		params.R_polaron_recombination = 1e10;
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
		// Check for heterojunction dependence using weakly donating/accepting bilayer
		params.R_polaron_recombination = 1e10;
		params.Enable_phase_restriction = false;
		params.Homo_acceptor = params.Homo_donor + 0.1;
		params.Lumo_acceptor = params.Lumo_acceptor + 0.1;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		double IQE6 = 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2.0 * (double)sim.getN_excitons_created());
		EXPECT_GT(IQE1, IQE6);
		// Check for higher order losses
		params.Height = 80;
		params.Thickness_donor = 40;
		params.Thickness_acceptor = 40;
		params.Triplet_lifetime_donor = 1e-7;
		params.Triplet_lifetime_acceptor = 1e-7;
		params.Enable_FRET_triplet_annihilation = true;
		params.R_exciton_exciton_annihilation_donor = 1e14;
		params.R_exciton_exciton_annihilation_acceptor = 1e14;
		params.R_exciton_polaron_annihilation_donor = 1e14;
		params.R_exciton_polaron_annihilation_acceptor = 1e14;
		params.R_exciton_isc_donor = 1e9;
		params.R_exciton_isc_acceptor = 1e9;
		params.R_exciton_risc_donor = 1e9;
		params.R_exciton_risc_acceptor = 1e9;
		params.E_exciton_ST_donor = 0.1;
		params.E_exciton_ST_acceptor = 0.1;
		params.R_polaron_recombination = 1e12;
		params.Enable_gaussian_polaron_delocalization = true;
		params.Polaron_delocalization_length = 4.0;
		params.Exciton_generation_rate_donor = 1e24;
		params.Exciton_generation_rate_donor = 1e24;
		params.Internal_potential = -1.0;
		params.N_tests = 1000;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		int N_exciton_exciton_annihilations1 = sim.getN_singlet_singlet_annihilations() + sim.getN_singlet_triplet_annihilations() + sim.getN_triplet_triplet_annihilations();
		int N_exciton_polaron_annihilations1 = sim.getN_singlet_polaron_annihilations() + sim.getN_triplet_polaron_annihilations();
		int N_bimolecular_recombinations1 = sim.getN_bimolecular_recombinations();
		// Check that higher order losses increase at higher generation rates
		params.Exciton_generation_rate_donor = 1e25;
		params.Exciton_generation_rate_donor = 1e25;
		sim = OSC_Sim();
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			success = sim.executeNextEvent();
			EXPECT_TRUE(success);
			if (!success) {
				cout << sim.getErrorMessage() << endl;
			}
		}
		int N_exciton_exciton_annihilations2 = sim.getN_singlet_singlet_annihilations() + sim.getN_singlet_triplet_annihilations() + sim.getN_triplet_triplet_annihilations();
		int N_exciton_polaron_annihilations2 = sim.getN_singlet_polaron_annihilations() + sim.getN_triplet_polaron_annihilations();
		int N_bimolecular_recombinations2 = sim.getN_bimolecular_recombinations();
		EXPECT_GT(N_exciton_exciton_annihilations2, N_exciton_exciton_annihilations1);
		EXPECT_GT(N_exciton_polaron_annihilations2, N_exciton_polaron_annihilations1);
		EXPECT_GT(N_bimolecular_recombinations2, N_bimolecular_recombinations1);
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
		double rate_constant = params.R_polaron_hopping_donor*exp(-2.0*params.Polaron_localization_donor);
		double expected_mobility = (rate_constant*1e-14) * (2.0 / 3.0) * (tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0)) * (1 / (K_b*params.Temperature));
		EXPECT_NEAR(expected_mobility, vector_avg(mobility_data), 1e-1*expected_mobility);
		// Hole ToF test with marcus hopping
		sim = OSC_Sim();
		params.Enable_miller_abrahams = false;
		params.Enable_marcus = true;
		EXPECT_TRUE(sim.init(params, 0));
		while (!sim.checkFinished()) {
			EXPECT_TRUE(sim.executeNextEvent());
		}
		transit_time_data = sim.getTransitTimeData();
		// Check that the transit time probability histogram sums to 1
		hist = sim.calculateTransitTimeHist(transit_time_data, (int)transit_time_data.size());
		cum_hist = calculateCumulativeHist(hist);
		EXPECT_NEAR(1.0, cum_hist.back().second, 1e-3);
		// Check the mobility compared to analytical expectation
		mobility_data = sim.calculateMobilityData(transit_time_data);
		rate_constant = (params.R_polaron_hopping_donor / sqrt(4.0*Pi*params.Reorganization_donor*K_b*params.Temperature))*exp(-2.0*params.Polaron_localization_donor)*exp(-intpow(params.Reorganization_donor + params.Internal_potential / params.Height, 2) / (4.0*params.Reorganization_donor*K_b*params.Temperature));
		expected_mobility = (rate_constant*1e-14) * (2.0 / 3.0) * (tgamma((dim + 1.0) / 2.0) / tgamma(dim / 2.0)) * (1 / (K_b*params.Temperature));
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
		params.Acceptor_conc = 0.99; // Dilute blend should behave like a neat electron transport material
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

	TEST_F(OSC_SimTest, InterfacialEnergyShiftTests) {
		// Test energy shift on bilayer without energetic disorder
		sim = OSC_Sim();
		auto params = params_default;
		params.Enable_neat = false;
		params.Enable_periodic_z = false;
		params.Enable_bilayer = true;
		params.Length = 100;
		params.Width = 100;
		params.Height = 40;
		params.Thickness_donor = 20;
		params.Thickness_acceptor = 20;
		params.Enable_interfacial_energy_shift = true;
		params.Energy_shift_donor = 0.01;
		params.Energy_shift_acceptor = 0.01;
		sim.init(params, 0);
		double expected_energy = params.Energy_shift_donor + (params.Energy_shift_donor * 4 / sqrt(2)) + (params.Energy_shift_donor * 4 / sqrt(3));
		EXPECT_DOUBLE_EQ(expected_energy, sim.getSiteEnergy(Coords(params.Length / 2, params.Width / 2, params.Height / 2)));
		// Test energy shift on bilayer with energetic disorder
		sim = OSC_Sim();
		params.Enable_gaussian_dos = true;
		params.Energy_stdev_donor = 0.05;
		params.Energy_stdev_acceptor = 0.05;
		sim.init(params, 0);
		vector<double> energies;
		for (int x = 0; x < params.Length; x++) {
			for (int y = 0; y < params.Width; y++) {
				energies.push_back(sim.getSiteEnergy(Coords(x, y, params.Height / 2)));
			}
		}
		double expected_energy_avg = params.Energy_shift_donor + (params.Energy_shift_donor * 4 / sqrt(2)) + (params.Energy_shift_donor * 4 / sqrt(3));
		EXPECT_NEAR(expected_energy_avg, vector_avg(energies), 5e-3);
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
