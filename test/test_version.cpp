// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "gtest/gtest.h"
#include "../src/Version.h"
#include "Utils.h"

using namespace std;
using namespace Utils;
using namespace Excimontec;

TEST(VersionTest, Version) {
	Version v1;
	v1.setVersion("v1.2-beta.4");
	Version v2;
	v2.setVersion("v1.2-beta.4");

	ASSERT_EQ(v1,v2);
	ASSERT_GE(v1,v2);
	ASSERT_LE(v1,v2);

	Version v3;
	v3.setVersion("v2.2-beta.4");
	ASSERT_NE(v1,v3);
	ASSERT_GE(v3,v1);
	ASSERT_GT(v3,v1);

	ASSERT_LE(v1,v3);
	ASSERT_LT(v1,v3);

	ASSERT_EQ(v1.getMajor(),1);
	ASSERT_EQ(v1.getMinor(),2);
	ASSERT_EQ(v1.getBugRelease(),4);
	string ver_type = v1.getReleaseType();
	ASSERT_STREQ("beta",ver_type.c_str());
	string ver = v1.getVersion();
	ASSERT_STREQ("v1.2-beta.4",ver.c_str());
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	// Redirect cout to NULL to suppress command line output during the tests
	cout.rdbuf(NULL);
	return RUN_ALL_TESTS();
}
