// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "Polaron.h"

using namespace std;

namespace Excimontec {

	// Initialize static class members
	const string Polaron::object_type = "Polaron";
	const string Polaron::Hop::event_type = "Polaron Hop";
	const string Polaron::Recombination::event_type = "Polaron Recombination";
	const string Polaron::Extraction::event_type = "Polaron Extraction";

}
