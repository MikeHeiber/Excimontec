// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "Exciton.h"

using namespace std;
using namespace Utils;
using namespace Excimontec;

// Initialize static class members
const string Exciton::object_type = "Exciton";
const string Exciton_Creation::event_type = "Exciton Creation";
const string Exciton_Hop::event_type = "Exciton Hop";
const string Exciton_Recombination::event_type = "Exciton Recombination";
const string Exciton_Dissociation::event_type = "Exciton Dissociation";
const string Exciton_Intersystem_Crossing::event_type = "Exciton Intersystem Crossing";
const string Exciton_Exciton_Annihilation::event_type = "Exciton-Exciton Annihilation";
const string Exciton_Polaron_Annihilation::event_type = "Exciton-Polaron Annihilation";

