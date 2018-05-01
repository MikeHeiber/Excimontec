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

namespace {

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
		mt19937 gen;
		vector<double> data;
		data.assign((int)3e7, 0);
		createGaussianDOSVector(data, 0.0, 0.15, gen);
		EXPECT_NEAR(0.0, vector_avg(data), 1e-4);
		EXPECT_NEAR(0.15, vector_stdev(data), 1e-4);
		auto hist = calculateProbabilityHist(data, 1000);
		EXPECT_NEAR(1.0, integrateData(hist), 1e-4);
		double peak = hist[499].second;
		EXPECT_NEAR(1.0 / sqrt(2.0*Pi*intpow(0.15, 2)), peak, 1e-2*peak);
		createGaussianDOSVector(data, 0.0, 0.05, gen);
		EXPECT_NEAR(0.0, vector_avg(data), 1e-4);
		EXPECT_NEAR(0.05, vector_stdev(data), 1e-4);
		hist = calculateProbabilityHist(data, 0.005);
		EXPECT_NEAR(1.0, integrateData(hist), 1e-4);
		vector<double> prob;
		for_each(hist.begin(), hist.end(), [&prob](pair<double, double>& x_y) {prob.push_back(x_y.second); });
		peak = *max_element(prob.begin(), prob.end());
		EXPECT_NEAR(1.0 / sqrt(2.0*Pi*intpow(0.05, 2)), peak, 1e-2*peak);
	}

	TEST(UtilsTests, ExponentialDOSTests) {
		mt19937 gen;
		vector<double> data;
		data.assign((int)2e7, 0);
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

namespace{
		
	class LatticeTest : public ::testing::Test {
		protected:
			mt19937 gen;
			Parameters_Lattice params_lattice;
			vector<Site> sites;
			Lattice lattice;

			void SetUp() {
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

namespace {

	class OSC_SimTest : public ::testing::Test {
		Parameters_OPV params_default;
		Parameters_Simulation params_base;
		OSC_Sim sim;
	};

}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
