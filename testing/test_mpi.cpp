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
			EXPECT_EQ(4.5, data_avg[0]);
			EXPECT_EQ(5.5, data_avg[1]);
			EXPECT_EQ(6.5, data_avg[2]);
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
