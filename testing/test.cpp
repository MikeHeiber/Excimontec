// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "gtest/gtest.h"
#include "OSC_Sim.h"
#include "Exciton.h"
#include "KMC_Lattice/Utils.h"

using namespace std;
using namespace Utils;

// Simple derived Simulation class
class TestSim : public Simulation {
public:
	bool checkFinished() const { return false; };
	bool executeNextEvent() { return false; };
};

namespace UtilsTests {

	TEST(UtilsTests, CoordsTests) {
		Coords coords{ 1, 1, 1 };
		EXPECT_EQ(1, coords.x);
		EXPECT_EQ(1, coords.y);
		EXPECT_EQ(1, coords.z);
		Coords coords2{ 0, 5, 10 };
		EXPECT_TRUE(coords != coords2);
		coords2.setXYZ(1, 1, 1);
		EXPECT_TRUE(coords == coords2);
	}

	TEST(UtilsTests, CalculateProbabilityHistTests) {
		vector<double> data{ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 };
		auto hist = calculateProbabilityHist(data, 1.0);
		EXPECT_EQ(3, (int)hist.size());
		EXPECT_DOUBLE_EQ((double)1 / (double)3, hist[0].second);
		EXPECT_DOUBLE_EQ((double)1 / (double)3, hist[1].second);
		EXPECT_DOUBLE_EQ((double)1 / (double)3, hist[2].second);
	}

	TEST(UtilsTests, ArrayStatsTests) {
		// positive ints
		int int_data[10];
		for (int i = 0; i < 10; i++) {
			int_data[i] = i + 1;
		}
		EXPECT_DOUBLE_EQ(5.5, array_avg(int_data, 10));
		EXPECT_NEAR(3.02765035409749, array_stdev(int_data, 10), 1e-14);
		// negative ints
		for (int i = 0; i < 10; i++) {
			int_data[i] = -(i + 1);
		}
		EXPECT_DOUBLE_EQ(-5.5, array_avg(int_data, 10));
		EXPECT_NEAR(3.02765035409749, array_stdev(int_data, 10), 1e-14);
		// positive doubles
		double double_data[10];
		for (int i = 0; i < 10; i++) {
			double_data[i] = 0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(2.75, array_avg(double_data, 10));
		EXPECT_NEAR(1.51382517704875, array_stdev(double_data, 10), 1e-14);
		// negative doubles
		for (int i = 0; i < 10; i++) {
			double_data[i] = -0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(-2.75, array_avg(double_data, 10));
		EXPECT_NEAR(1.51382517704875, array_stdev(double_data, 10), 1e-14);
	}
	TEST(UtilsTests, VectorStatsTests) {
		// positive ints
		vector<int> int_data;
		int_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			int_data[i] = i + 1;
		}
		EXPECT_DOUBLE_EQ(5.5, vector_avg(int_data));
		EXPECT_NEAR(3.02765035409749, vector_stdev(int_data), 1e-14);
		// negative ints
		for (int i = 0; i < 10; i++) {
			int_data[i] = -(i + 1);
		}
		EXPECT_DOUBLE_EQ(-5.5, vector_avg(int_data));
		EXPECT_NEAR(3.02765035409749, vector_stdev(int_data), 1e-14);
		// positive doubles
		vector<double> double_data;
		double_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			double_data[i] = 0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(2.75, vector_avg(double_data));
		EXPECT_NEAR(1.51382517704875, vector_stdev(double_data), 1e-14);
		// negative doubles
		double_data.assign(10, 0);
		for (int i = 0; i < 10; i++) {
			double_data[i] = -0.5*(i + 1);
		}
		EXPECT_DOUBLE_EQ(-2.75, vector_avg(double_data));
		EXPECT_NEAR(1.51382517704875, vector_stdev(double_data), 1e-14);
	}

	TEST(UtilsTests, IntPowTests) {
		EXPECT_DOUBLE_EQ(1.0, intpow(2.5, 0));
		EXPECT_DOUBLE_EQ(2.5, intpow(2.5, 1));
		EXPECT_DOUBLE_EQ(6.25, intpow(2.5, 2));
		EXPECT_DOUBLE_EQ(9536.7431640625, intpow(2.5, 10));
		EXPECT_DOUBLE_EQ(0.4, intpow(2.5, -1));
		EXPECT_DOUBLE_EQ(0.16, intpow(2.5, -2));
		EXPECT_DOUBLE_EQ(1.048576e-4, intpow(2.5, -10));
		EXPECT_DOUBLE_EQ(1.0, intpow(15.04564, 0));
		EXPECT_DOUBLE_EQ(1e-21, intpow(1e-7, 3));
	}

	TEST(UtilsTests, GaussianDOSTests) {
		mt19937 gen(std::random_device{}());
		vector<double> data((int)2e7,0.0);
		createGaussianDOSVector(data, 0.0, 0.15, gen);
		EXPECT_NEAR(0.0, vector_avg(data), 1e-4);
		EXPECT_NEAR(0.15, vector_stdev(data), 1e-4);
		auto hist = calculateProbabilityHist(data, 1000);
		EXPECT_NEAR(1.0, integrateData(hist), 1e-4);
		double peak = hist[499].second;
		EXPECT_NEAR(1.0 / sqrt(2.0*Pi*intpow(0.15, 2)), peak, 5e-2*peak);
		createGaussianDOSVector(data, 0.0, 0.05, gen);
		EXPECT_NEAR(0.0, vector_avg(data), 1e-4);
		EXPECT_NEAR(0.05, vector_stdev(data), 1e-4);
		hist = calculateProbabilityHist(data, 0.005);
		EXPECT_NEAR(1.0, integrateData(hist), 1e-4);
		vector<double> prob;
		for_each(hist.begin(), hist.end(), [&prob](pair<double, double>& x_y) {prob.push_back(x_y.second); });
		peak = *max_element(prob.begin(), prob.end());
		EXPECT_NEAR(1.0 / sqrt(2.0*Pi*intpow(0.05, 2)), peak, 5e-2*peak);
	}

	TEST(UtilsTests, ExponentialDOSTests) {
		mt19937 gen(std::random_device{}());
		vector<double> data((int)2e7, 0.0);
		createExponentialDOSVector(data, 0.0, 0.1, gen);
		auto hist = calculateProbabilityHist(data, 1000);
		EXPECT_NEAR(1.0, integrateData(hist), 1e-4);
		vector<double> prob;
		for_each(hist.begin(), hist.end(), [&prob](pair<double, double>& x_y) {prob.push_back(x_y.second); });
		double peak = *max_element(prob.begin(), prob.end());
		EXPECT_NEAR(0.5*(1.0 / 0.1), peak, 1e-2*peak);
		createExponentialDOSVector(data, 0.0, 0.05, gen);
		hist = calculateProbabilityHist(data, 1000);
		EXPECT_NEAR(1.0, integrateData(hist), 1e-4);
		prob.clear();
		for_each(hist.begin(), hist.end(), [&prob](pair<double, double>& x_y) {prob.push_back(x_y.second); });
		peak = *max_element(prob.begin(), prob.end());
		EXPECT_NEAR(0.5*(1.0 / 0.05), peak, 1e-2*peak);
	}

	TEST(UtilsTests, ImportBooleanTests) {
		bool error_status;
		EXPECT_TRUE(importBooleanParam("true", error_status));
		EXPECT_TRUE(importBooleanParam(" true  ", error_status));
		EXPECT_FALSE(importBooleanParam("false", error_status));
		EXPECT_FALSE(importBooleanParam("   false", error_status));
		EXPECT_FALSE(importBooleanParam("blah", error_status));
		EXPECT_FALSE(importBooleanParam("blah  ", error_status));
		EXPECT_TRUE(error_status);
	}

	TEST(UtilsTests, RemoveDuplicatesTests) {
		vector<int> vec{ 0, 1, 1, 2, 3, 1, 4, 2 };
		removeDuplicates(vec);
		EXPECT_EQ(2, vec[2]);
		EXPECT_EQ(5, (int)vec.size());
		vec = { 0, 1, 2, 1 };
		removeDuplicates(vec);
		EXPECT_EQ(2, vec[2]);
		EXPECT_EQ(3, (int)vec.size());
		vector<double> vec2{ 0.0, 1.0, 1.0, 2.0, 3.0, 1.0, 4.0, 2.0 };
		removeDuplicates(vec2);
		EXPECT_DOUBLE_EQ(2.0, vec2[2]);
		EXPECT_EQ(5, (int)vec2.size());
		Coords coords1{ 1,2,3 };
		Coords coords2{ 1,2,3 };
		Coords coords3{ 4,5,6 };
		vector<Coords> vec3{ coords1, coords1, coords2, coords3, coords1, coords3, coords2 };
		removeDuplicates(vec3);
		EXPECT_EQ(4, vec3[1].x);
		EXPECT_EQ(2, (int)vec3.size());
		Exciton exciton1(0, 1, { 0,0,0 });
		Exciton exciton2(0, 2, { 1,1,1 });
		vector<Exciton*> object_ptrs{ &exciton1, &exciton1, &exciton2, &exciton1, &exciton2 };
		removeDuplicates(object_ptrs);
		EXPECT_EQ(2, (int)object_ptrs.size());
		vector<int> vec4 = {};
		removeDuplicates(vec4);
		EXPECT_EQ(0, (int)vec4.size());
		vector<int> vec5 = { 0 };
		removeDuplicates(vec5);
		EXPECT_EQ(1, (int)vec5.size());
		vector<int> vec6 = { 0,0 };
		removeDuplicates(vec6);
		EXPECT_EQ(1, (int)vec6.size());
	}

}

namespace LatticeTests{
		
	class LatticeTest : public ::testing::Test {
		protected:
			mt19937 gen;
			Parameters_Lattice params_lattice;
			vector<Site> sites;
			Lattice lattice;

			void SetUp() {
				gen.seed(std::random_device{}());
				// Setup params
				params_lattice.Enable_periodic_x = true;
				params_lattice.Enable_periodic_y = true;
				params_lattice.Enable_periodic_z = true;
				params_lattice.Length = 50;
				params_lattice.Width = 50;
				params_lattice.Height = 50;
				params_lattice.Unit_size = 1.0;
				// Initialize Lattice object
				lattice.init(params_lattice, &gen);
				Site site;
				sites.assign(lattice.getNumSites(), site);
				vector<Site*> site_ptrs((int)sites.size());
				for (int i = 0; i < (int)sites.size(); i++) {
					site_ptrs[i] = &sites[i];
				}
				lattice.setSitePointers(site_ptrs);
			}
	};

	TEST_F(LatticeTest, InitializationTests) {
		EXPECT_EQ(50, lattice.getLength());
		EXPECT_EQ(50, lattice.getWidth());
		EXPECT_EQ(50, lattice.getHeight());
		EXPECT_DOUBLE_EQ(1.0, lattice.getUnitSize());
		EXPECT_EQ((long int)50 * 50 * 50, lattice.getNumSites());
	}

	TEST_F(LatticeTest, CalculateDestCoordsTests) {
		Coords coords_i{ 49, 49, 49 };
		Coords coords_f;
		lattice.calculateDestinationCoords(coords_i, 2, -2, 0, coords_f);
		EXPECT_EQ(1, coords_f.x);
		EXPECT_EQ(47, coords_f.y);
		EXPECT_EQ(49, coords_f.z);
		coords_i.setXYZ(0, 0, 49);
		lattice.calculateDestinationCoords(coords_i, 1, -2, 2, coords_f);
		EXPECT_EQ(1, coords_f.x);
		EXPECT_EQ(48, coords_f.y);
		EXPECT_EQ(1, coords_f.z);
	}

	TEST_F(LatticeTest, PeriodicCrossingTests) {
		Coords coords_i{ 49, 49, 49 };
		Coords coords_f{ 0, 49, 48 };
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), -50);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), 0);
		coords_i.setXYZ(0, 0, 49);
		coords_f.setXYZ(1, 49, 0);
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 50);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), -50);
		coords_i.setXYZ(4, 5, 6);
		coords_f.setXYZ(3, 6, 5);
		EXPECT_EQ(lattice.calculateDX(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDY(coords_i, coords_f), 0);
		EXPECT_EQ(lattice.calculateDZ(coords_i, coords_f), 0);
	}

	TEST_F(LatticeTest, CheckMoveValidityTests) {
		params_lattice.Enable_periodic_x = false;
		params_lattice.Enable_periodic_y = false;
		params_lattice.Enable_periodic_z = true;
		Lattice lattice2;
		lattice2.init(params_lattice, &gen);
		Coords coords1{ 49, 0, 49 };
		EXPECT_FALSE(lattice2.checkMoveValidity(coords1, 1, 0, 0));
		EXPECT_FALSE(lattice2.checkMoveValidity(coords1, 0, -1, 0));
		EXPECT_TRUE(lattice2.checkMoveValidity(coords1, 0, 0, 1));
		EXPECT_TRUE(lattice2.checkMoveValidity(coords1, -1, 1, -1));
	}

	TEST_F(LatticeTest, RandomSiteGenTests) {
		Coords coords;
		int N_test = 40000000;
		vector<int> xcoords(N_test);
		vector<int> ycoords(N_test);
		vector<int> zcoords(N_test);
		for (int i = 0; i < N_test; i++) {
			coords = lattice.generateRandomCoords();
			EXPECT_TRUE(coords.x >= 0 && coords.x < lattice.getLength());
			EXPECT_TRUE(coords.y >= 0 && coords.y < lattice.getWidth());
			EXPECT_TRUE(coords.z >= 0 && coords.z < lattice.getHeight());
			xcoords[i] = coords.x;
			ycoords[i] = coords.y;
			zcoords[i] = coords.z;
		}
		EXPECT_NEAR(24.5, vector_avg(xcoords), 2e-2);
		EXPECT_NEAR(lattice.getLength() / sqrt(12.0), vector_stdev(xcoords), 1e-2);
		EXPECT_NEAR(24.5, vector_avg(ycoords), 2e-2);
		EXPECT_NEAR(lattice.getWidth() / sqrt(12.0), vector_stdev(ycoords), 1e-2);
		EXPECT_NEAR(24.5, vector_avg(zcoords), 2e-2);
		EXPECT_NEAR(lattice.getHeight() / sqrt(12.0), vector_stdev(zcoords), 1e-2);
	}

	TEST_F(LatticeTest, LatticeDistanceTests) {
		Coords coords_i{ 49, 49, 49 };
		Coords coords_f{ 1, 49, 47 };
		EXPECT_EQ(8, lattice.calculateLatticeDistanceSquared(coords_i, coords_f));
		coords_i.setXYZ(4, 5, 6);
		coords_f.setXYZ(3, 6, 5);
		EXPECT_EQ(3, lattice.calculateLatticeDistanceSquared(coords_i, coords_f));
		coords_i.setXYZ(0, 49, 1);
		coords_f.setXYZ(48, 1, 49);
		EXPECT_EQ(12, lattice.calculateLatticeDistanceSquared(coords_i, coords_f));
	}

	TEST_F(LatticeTest, GetSiteCoordsTests) {
		Coords coords1, coords2;
		int index;
		for (int i = 0; i < 100; i++) {
			coords1 = lattice.generateRandomCoords();
			index = lattice.getSiteIndex(coords1);
			coords2 = lattice.getSiteCoords(index);
			EXPECT_EQ(coords1.x, coords2.x);
			EXPECT_EQ(coords1.y, coords2.y);
			EXPECT_EQ(coords1.z, coords2.z);
		}
	}

	TEST_F(LatticeTest, OccupancyTests) {
		Coords coords;
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					coords.setXYZ(x, y, z);
					EXPECT_FALSE(lattice.isOccupied(coords));
				}
			}
		}
		coords.setXYZ(5, 5, 5);
		lattice.setOccupied(coords);
		EXPECT_TRUE(lattice.isOccupied(coords));
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					if (lattice.isOccupied(coords)) {
						EXPECT_EQ(5, coords.x);
						EXPECT_EQ(5, coords.y);
						EXPECT_EQ(5, coords.z);
					}
				}
			}
		}
		lattice.clearOccupancy(coords);
		EXPECT_FALSE(lattice.isOccupied(coords));
	}

	TEST_F(LatticeTest, SetSitePointersTests) {
		Site site;
		vector<Site*> site_ptrs(10, &site);
		EXPECT_FALSE(lattice.setSitePointers(site_ptrs));
		site_ptrs.assign(50 * 50 * 50, &site);
		EXPECT_TRUE(lattice.setSitePointers(site_ptrs));
	}
}

namespace EventTests {

	class EventTest : public ::testing::Test {
	protected:
		Parameters_OPV params_default;
		Parameters_Simulation params_base;
		TestSim test_sim;
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
			params_default.Length = 50;
			params_default.Width = 50;
			params_default.Height = 50;
			params_default.Unit_size = 1.0;
			params_default.Temperature = 300;
			// Output files
			params_default.Logfile = NULL;
			// Additional General Parameters
			params_default.Internal_potential = 0.0;
			// Morphology Parameters
			params_default.Enable_neat = true;
			params_default.Enable_bilayer = false;
			params_default.Thickness_donor = 25;
			params_default.Thickness_acceptor = 25;
			params_default.Enable_random_blend = false;
			params_default.Acceptor_conc = 0.5;
			params_default.Enable_import_morphology = false;
			params_default.Morphology_file = NULL;
			// Test Parameters
			params_default.N_tests = 100;
			params_default.Enable_exciton_diffusion_test = true;
			params_default.Enable_ToF_test = false;
			params_default.ToF_polaron_type = false;
			params_default.ToF_initial_polarons = 1;
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
			params_default.Exciton_generation_rate_donor = 1e22;
			params_default.Exciton_generation_rate_acceptor = 1e22;
			params_default.Singlet_lifetime_donor = 500e-12;
			params_default.Singlet_lifetime_acceptor = 500e-12;
			params_default.Triplet_lifetime_donor = 1e-6;
			params_default.Triplet_lifetime_acceptor = 1e-6;
			params_default.R_singlet_hopping_donor = 1e12;
			params_default.R_singlet_hopping_acceptor = 1e12;
			params_default.Singlet_localization_donor = 1.0;
			params_default.Singlet_localization_acceptor = 1.0;
			params_default.R_triplet_hopping_donor = 1e12;
			params_default.R_triplet_hopping_acceptor = 1e12;
			params_default.Triplet_localization_donor = 2.0;
			params_default.Triplet_localization_acceptor = 2.0;
			params_default.Enable_FRET_triplet_annihilation = false;
			params_default.R_exciton_exciton_annihilation_donor = 1e12;
			params_default.R_exciton_exciton_annihilation_acceptor = 1e12;
			params_default.R_exciton_polaron_annihilation_donor = 1e12;
			params_default.R_exciton_polaron_annihilation_acceptor = 1e12;
			params_default.FRET_cutoff = 3;
			params_default.E_exciton_binding_donor = 0.5;
			params_default.E_exciton_binding_acceptor = 0.5;
			params_default.R_exciton_dissociation_donor = 1e14;
			params_default.R_exciton_dissociation_acceptor = 1e14;
			params_default.Exciton_dissociation_cutoff = 3;
			params_default.R_exciton_isc_donor = 1e7;
			params_default.R_exciton_isc_acceptor = 1e7;
			params_default.R_exciton_risc_donor = 1e7;
			params_default.R_exciton_risc_acceptor = 1e7;
			params_default.E_exciton_ST_donor = 0.7;
			params_default.E_exciton_ST_acceptor = 0.7;
			// Polaron Parameters
			params_default.Enable_phase_restriction = true;
			params_default.R_polaron_hopping_donor = 1e12;
			params_default.R_polaron_hopping_acceptor = 1e12;
			params_default.Polaron_localization_donor = 2.0;
			params_default.Polaron_localization_acceptor = 2.0;
			params_default.Enable_miller_abrahams = true;
			params_default.Enable_marcus = false;
			params_default.Reorganization_donor = 0.2;
			params_default.Reorganization_acceptor = 0.2;
			params_default.R_polaron_recombination = 1e12;
			params_default.Polaron_hopping_cutoff = 3;
			params_default.Enable_gaussian_polaron_delocalization = false;
			params_default.Polaron_delocalization_length = 1.0;
			// Additional Lattice Parameters
			params_default.Homo_donor = 5.0;
			params_default.Lumo_donor = 3.0;
			params_default.Homo_acceptor = 6.0;
			params_default.Lumo_acceptor = 4.0;
			params_default.Enable_gaussian_dos = true;
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
			params_default.Coulomb_cutoff = 15;
			// Initialize TestSim object
			params_base = params_default;
			test_sim.init(params_base, 0);
		}
	};

	TEST_F(EventTest, CalculateExecutionTimeTests) {
		Event event(&test_sim);
		// Generate collection of wait times
		vector<double> times((int)3e7);
		for (int i = 0; i < (int)times.size(); i++) {
			event.calculateExecutionTime(1e9);
			times[i] = event.getExecutionTime();
		}
		// Calculate probability distribution of calculated times
		auto hist = calculateProbabilityHist(times, 1000);
		// Test probability distribution
		EXPECT_NEAR(1.0, integrateData(hist), 2e-2);
		EXPECT_NEAR(1e9, hist[0].second, 2e7);
		auto it = find_if(hist.begin(), hist.end(), [](pair<double, double>& x_y) {return x_y.first > 1e-9; });
		it--;
		EXPECT_NEAR(1e9 / 2.7182818284590, (*it).second, 2e7);
		// Test average wait time
		EXPECT_NEAR(vector_avg(times), 1e-9, 2e-12);
		// Generate collection of wait times
		for (int i = 0; i < (int)times.size(); i++) {
			event.calculateExecutionTime(1e12);
			times[i] = event.getExecutionTime();
		}
		// Calculate probability distribution of calculated times
		hist = calculateProbabilityHist(times, 1000);
		// Test probability distribution
		EXPECT_NEAR(1.0, integrateData(hist), 2e-2);
		EXPECT_NEAR(1e12, hist[0].second, 2e10);
		it = find_if(hist.begin(), hist.end(), [](pair<double, double>& x_y) {return x_y.first > 1e-12; });
		it--;
		EXPECT_NEAR(1e12 / 2.7182818284590, (*it).second, 2e10);
		// Test average wait time
		EXPECT_NEAR(vector_avg(times), 1e-12, 2e-15);
	}
}

namespace ObjectTests {
	class EventTest : public ::testing::Test {
	protected:
		Parameters_OPV params_default;
		Parameters_Simulation params_base;
		TestSim test_sim;
		void SetUp() {
			// Initialize default parameters
				// General Parameters
			params_default.Enable_FRM = false;
			params_default.Enable_selective_recalc = true;
			params_default.Recalc_cutoff = 3;
			params_default.Enable_full_recalc = false;
			params_default.Enable_logging = false;
			params_default.Enable_periodic_x = true;
			params_default.Enable_periodic_y = true;
			params_default.Enable_periodic_z = true;
			params_default.Length = 50;
			params_default.Width = 50;
			params_default.Height = 50;
			params_default.Unit_size = 1.0;
			params_default.Temperature = 300;
			// Output files
			params_default.Logfile = NULL;
			// Additional General Parameters
			params_default.Internal_potential = 0.0;
			// Morphology Parameters
			params_default.Enable_neat = true;
			params_default.Enable_bilayer = false;
			params_default.Thickness_donor = 25;
			params_default.Thickness_acceptor = 25;
			params_default.Enable_random_blend = false;
			params_default.Acceptor_conc = 0.5;
			params_default.Enable_import_morphology = false;
			params_default.Morphology_file = NULL;
			// Test Parameters
			params_default.N_tests = 100;
			params_default.Enable_exciton_diffusion_test = true;
			params_default.Enable_ToF_test = false;
			params_default.ToF_polaron_type = false;
			params_default.ToF_initial_polarons = 1;
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
			params_default.Exciton_generation_rate_donor = 1e22;
			params_default.Exciton_generation_rate_acceptor = 1e22;
			params_default.Singlet_lifetime_donor = 500e-12;
			params_default.Singlet_lifetime_acceptor = 500e-12;
			params_default.Triplet_lifetime_donor = 1e-6;
			params_default.Triplet_lifetime_acceptor = 1e-6;
			params_default.R_singlet_hopping_donor = 1e12;
			params_default.R_singlet_hopping_acceptor = 1e12;
			params_default.Singlet_localization_donor = 1.0;
			params_default.Singlet_localization_acceptor = 1.0;
			params_default.R_triplet_hopping_donor = 1e12;
			params_default.R_triplet_hopping_acceptor = 1e12;
			params_default.Triplet_localization_donor = 2.0;
			params_default.Triplet_localization_acceptor = 2.0;
			params_default.Enable_FRET_triplet_annihilation = false;
			params_default.R_exciton_exciton_annihilation_donor = 1e12;
			params_default.R_exciton_exciton_annihilation_acceptor = 1e12;
			params_default.R_exciton_polaron_annihilation_donor = 1e12;
			params_default.R_exciton_polaron_annihilation_acceptor = 1e12;
			params_default.FRET_cutoff = 3;
			params_default.E_exciton_binding_donor = 0.5;
			params_default.E_exciton_binding_acceptor = 0.5;
			params_default.R_exciton_dissociation_donor = 1e14;
			params_default.R_exciton_dissociation_acceptor = 1e14;
			params_default.Exciton_dissociation_cutoff = 3;
			params_default.R_exciton_isc_donor = 1e7;
			params_default.R_exciton_isc_acceptor = 1e7;
			params_default.R_exciton_risc_donor = 1e7;
			params_default.R_exciton_risc_acceptor = 1e7;
			params_default.E_exciton_ST_donor = 0.7;
			params_default.E_exciton_ST_acceptor = 0.7;
			// Polaron Parameters
			params_default.Enable_phase_restriction = true;
			params_default.R_polaron_hopping_donor = 1e12;
			params_default.R_polaron_hopping_acceptor = 1e12;
			params_default.Polaron_localization_donor = 2.0;
			params_default.Polaron_localization_acceptor = 2.0;
			params_default.Enable_miller_abrahams = true;
			params_default.Enable_marcus = false;
			params_default.Reorganization_donor = 0.2;
			params_default.Reorganization_acceptor = 0.2;
			params_default.R_polaron_recombination = 1e12;
			params_default.Polaron_hopping_cutoff = 3;
			params_default.Enable_gaussian_polaron_delocalization = false;
			params_default.Polaron_delocalization_length = 1.0;
			// Additional Lattice Parameters
			params_default.Homo_donor = 5.0;
			params_default.Lumo_donor = 3.0;
			params_default.Homo_acceptor = 6.0;
			params_default.Lumo_acceptor = 4.0;
			params_default.Enable_gaussian_dos = true;
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
			params_default.Coulomb_cutoff = 15;
			// Initialize TestSim object
			params_base = params_default;
			test_sim.init(params_base, 0);
		}
	};

	TEST(ObjectTests, CalculateDisplacementTests) {
		// Assume object is created in a 50 x 50 x 50 lattice
		Object object1(0, 1, { 0,0,0 });
		object1.setCoords({ 0,0,5 });
		EXPECT_DOUBLE_EQ(5.0, object1.calculateDisplacement());
		object1.setCoords({ 0,0,49 });
		object1.incrementDZ(-50);
		EXPECT_DOUBLE_EQ(1.0, object1.calculateDisplacement());
		object1.setCoords({ 49,49,49 });
		object1.incrementDX(-50);
		object1.incrementDY(-50);
		EXPECT_DOUBLE_EQ(sqrt(3), object1.calculateDisplacement());
		object1.setCoords({ 49,0,49 });
		object1.incrementDY(50);
		EXPECT_DOUBLE_EQ(sqrt(2), object1.calculateDisplacement());
	}
}

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
			params_default.Length = 50;
			params_default.Width = 50;
			params_default.Height = 50;
			params_default.Unit_size = 1.0;
			params_default.Temperature = 300;
			// Output files
			params_default.Logfile = NULL;
			// Additional General Parameters
			params_default.Internal_potential = 0.0;
			// Morphology Parameters
			params_default.Enable_neat = true;
			params_default.Enable_bilayer = false;
			params_default.Thickness_donor = 25;
			params_default.Thickness_acceptor = 25;
			params_default.Enable_random_blend = false;
			params_default.Acceptor_conc = 0.5;
			params_default.Enable_import_morphology = false;
			params_default.Morphology_file = NULL;
			// Test Parameters
			params_default.N_tests = 100;
			params_default.Enable_exciton_diffusion_test = true;
			params_default.Enable_ToF_test = false;
			params_default.ToF_polaron_type = false;
			params_default.ToF_initial_polarons = 1;
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
			params_default.Exciton_generation_rate_donor = 1e22;
			params_default.Exciton_generation_rate_acceptor = 1e22;
			params_default.Singlet_lifetime_donor = 500e-12;
			params_default.Singlet_lifetime_acceptor = 500e-12;
			params_default.Triplet_lifetime_donor = 1e-6;
			params_default.Triplet_lifetime_acceptor = 1e-6;
			params_default.R_singlet_hopping_donor = 1e12;
			params_default.R_singlet_hopping_acceptor = 1e12;
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
			params_default.FRET_cutoff = 3;
			params_default.E_exciton_binding_donor = 0.5;
			params_default.E_exciton_binding_acceptor = 0.5;
			params_default.R_exciton_dissociation_donor = 1e14;
			params_default.R_exciton_dissociation_acceptor = 1e14;
			params_default.Exciton_dissociation_cutoff = 3;
			params_default.R_exciton_isc_donor = 1e7;
			params_default.R_exciton_isc_acceptor = 1e7;
			params_default.R_exciton_risc_donor = 1e7;
			params_default.R_exciton_risc_acceptor = 1e7;
			params_default.E_exciton_ST_donor = 0.7;
			params_default.E_exciton_ST_acceptor = 0.7;
			// Polaron Parameters
			params_default.Enable_phase_restriction = true;
			params_default.R_polaron_hopping_donor = 1e12;
			params_default.R_polaron_hopping_acceptor = 1e12;
			params_default.Polaron_localization_donor = 2.0;
			params_default.Polaron_localization_acceptor = 2.0;
			params_default.Enable_miller_abrahams = true;
			params_default.Enable_marcus = false;
			params_default.Reorganization_donor = 0.2;
			params_default.Reorganization_acceptor = 0.2;
			params_default.R_polaron_recombination = 1e12;
			params_default.Polaron_hopping_cutoff = 3;
			params_default.Enable_gaussian_polaron_delocalization = false;
			params_default.Polaron_delocalization_length = 1.0;
			// Additional Lattice Parameters
			params_default.Homo_donor = 5.0;
			params_default.Lumo_donor = 3.0;
			params_default.Homo_acceptor = 6.0;
			params_default.Lumo_acceptor = 4.0;
			params_default.Enable_gaussian_dos = true;
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
			params_default.Coulomb_cutoff = 15;

			// Initialize OSC_Sim object
			sim.init(params_default, 0);
		}
	};

	TEST_F(OSC_SimTest, ParameterCheckTests) {
		Parameters_OPV params = params_default;
		EXPECT_TRUE(sim.init(params, 0));
		params.Enable_bilayer = true;
		EXPECT_FALSE(sim.init(params, 0));
		params.Enable_bilayer = false;
		params.Enable_ToF_test = true;
		EXPECT_FALSE(sim.init(params, 0));
	}

	TEST_F(OSC_SimTest, CorrelatedDisorderGaussianKernelTests) {
		Parameters_OPV params = params_default;
		params.Energy_stdev_donor = 0.05;
		params.Energy_stdev_acceptor = 0.05;
		params.Enable_correlated_disorder = true;
		params.Enable_gaussian_kernel = true;
		params.Enable_power_kernel = false;
		// correlation length = 1.2, Unit_size = 1.0
		params.Disorder_correlation_length = 1.2;
		sim.init(params, 0);
		auto energies = sim.getSiteEnergies(1);
		EXPECT_NEAR(0.0, vector_avg(energies), 5e-3);
		EXPECT_NEAR(0.05, vector_stdev(energies), 1e-3);
		auto correlation_data = sim.getDOSCorrelationData();
		EXPECT_NEAR(1/exp(1), interpolateData(correlation_data, 1.2), 0.05);
		// correlation length = 1.1, Unit_size = 0.8
		params.Unit_size = 0.8;
		params.Length = 40;
		params.Width = 40;
		params.Height = 40;
		params.Disorder_correlation_length = 1.1;
		sim.init(params, 0);
		energies = sim.getSiteEnergies(1);
		EXPECT_NEAR(0.0, vector_avg(energies), 5e-3);
		EXPECT_NEAR(0.05, vector_stdev(energies), 1e-3);
		correlation_data = sim.getDOSCorrelationData();
		EXPECT_NEAR(1/exp(1), interpolateData(correlation_data, 1.1), 0.05);
		// correlation length = 1.4, Unit_size = 1.2
		params.Unit_size = 1.2;
		params.Length = 60;
		params.Width = 60;
		params.Height = 60;
		params.Disorder_correlation_length = 1.4;
		sim.init(params, 0);
		energies = sim.getSiteEnergies(1);
		EXPECT_NEAR(0.0, vector_avg(energies), 5e-3);
		EXPECT_NEAR(0.05, vector_stdev(energies), 1e-3);
		correlation_data = sim.getDOSCorrelationData();
		EXPECT_NEAR(1 / exp(1), interpolateData(correlation_data, 1.4), 0.05);
	}

}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
