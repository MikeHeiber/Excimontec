// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "gtest/gtest.h"
#include "KMC_Lattice/Utils.h"
#include <mpi.h>

using namespace std;
using namespace Utils;

namespace MPI_Tests {

	class MPI_Test : public ::testing::Test {
	protected:
		int procid;
		int nproc;
		void SetUp() {
			MPI_Comm_size(MPI_COMM_WORLD, &nproc);
			MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		}
	};

	TEST_F(MPI_Test, CalculateVectorAvgTests) {
		// Seteup unique data vectors on each proc
		vector<double> data(3);
		for (int i = 0; i < (int)data.size(); i++) {
			data[i] = 3 * procid + i;
		}
		vector<double> data_avg = MPI_calculateVectorAvg(data);
		if (procid == 0) {
			EXPECT_EQ(3, (int)data_avg.size());
			EXPECT_DOUBLE_EQ(4.5, data_avg[0]);
			EXPECT_DOUBLE_EQ(5.5, data_avg[1]);
			EXPECT_DOUBLE_EQ(6.5, data_avg[2]);
		}
	}

	TEST_F(MPI_Test, CalculateVectorSumTests) {
		// Seteup unique data vectors on each proc
		vector<double> data(3);
		for (int i = 0; i < (int)data.size(); i++) {
			data[i] = 3 * procid + i;
		}
		vector<double> data_sum = MPI_calculateVectorSum(data);
		if (procid == 0) {
			EXPECT_EQ(3, (int)data_sum.size());
			EXPECT_DOUBLE_EQ(18.0, data_sum[0]);
			EXPECT_DOUBLE_EQ(22.0, data_sum[1]);
			EXPECT_DOUBLE_EQ(26.0, data_sum[2]);
		}
		vector<int> data2(3);
		for (int i = 0; i < (int)data2.size(); i++) {
			data2[i] = 3 * procid + i;
		}
		vector<int> data_sum2 = MPI_calculateVectorSum(data2);
		if (procid == 0) {
			EXPECT_EQ(3, (int)data_sum2.size());
			EXPECT_EQ(18, data_sum2[0]);
			EXPECT_EQ(22, data_sum2[1]);
			EXPECT_EQ(26, data_sum2[2]);
		}
	}

	TEST_F(MPI_Test, GatherVectorsTests) {
		// Seteup unique data vectors on each proc
		vector<double> data(3);
		for (int i = 0; i < (int)data.size(); i++) {
			data[i] = 3 * procid + i;
		}
		vector<double> data_all = MPI_gatherVectors(data);
		if (procid == 0) {
			EXPECT_EQ(12, (int)data_all.size());
			for (int i = 0; i < 3 * nproc; i++) {
				EXPECT_DOUBLE_EQ((double)i, data_all[i]);
			}
		}
		vector<int> data2(3);
		for (int i = 0; i < (int)data2.size(); i++) {
			data2[i] = 3 * procid + i;
		}
		vector<int> data_all2 = MPI_gatherVectors(data2);
		if (procid == 0) {
			EXPECT_EQ(12, (int)data_all2.size());
			for (int i = 0; i < 3 * nproc; i++) {
				EXPECT_EQ(i, data_all2[i]);
			}
		}
	}
}

int main(int argc, char **argv) {
	int result = 0;
	::testing::InitGoogleTest(&argc, argv);
	MPI_Init(&argc, &argv);
	result = RUN_ALL_TESTS();
	MPI_Finalize();
	return result;
}
