// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "gtest/gtest.h"
#include "Version.h"
#include "Utils.h"

using namespace std;
using namespace Utils;
using namespace Excimontec;

namespace Version {

	TEST_F(VersionTest, Version) {
		Version v1("v1.2-beta.4");
		Version v2("v1.2-beta.4");

		ASSERT_EQ(v1,v2);
		ASSERT_Ge(v1,v2);
		ASSERT_Le(v1,v2);

		Version v3("v2.2-beta.4");
		ASSERT_Ne(v1,v3);
		ASSERT_Ge(v3,v1);
		ASSERT_Gt(v3,v1);
		
		ASSERT_Le(v1,v3);
		ASSERT_Lt(v1,v3);

		ASSERT_EQ(v1.getMajor(),1);
		ASSERT_EQ(v1.getMinor(),2);
		ASSERT_EQ(v1.getBugRelease(),4);
		ASSERT_STEQ(v1.getReleaseType(),"beta");
		ASSERT_STEQ(v1.getVersion(),"v1.2-beta.4");
	
	}
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	// Redirect cout to NULL to suppress command line output during the tests
	cout.rdbuf(NULL);
	return RUN_ALL_TESTS();
}
