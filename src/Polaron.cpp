// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "Polaron.h"

using namespace std;
using namespace Utils;
using namespace Excimontec;

// Initialize static class members
const string Excimontec::Polaron::object_type = "Polaron";
const string Excimontec::Polaron_Hop::event_type = "Polaron Hop";
const string Excimontec::Polaron_Recombination::event_type = "Polaron Recombination";
const string Excimontec::Polaron_Extraction::event_type = "Polaron Extraction";
