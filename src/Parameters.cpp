// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "Parameters.h"

using namespace std;
using namespace KMC_Lattice;

namespace Excimontec {

	bool Parameters::checkParameters() const {
		// Check lattice parameters and other general parameters
		if(!Parameters_Simulation::checkParameters()) {
			return false;
		}
		if (Enable_selective_recalc && Recalc_cutoff < FRET_cutoff) {
			cout << "Error! The event recalculation cutoff radius must not be less than the FRET cutoff radius." << endl;
			return false;
		}
		if (Enable_selective_recalc && Recalc_cutoff < Polaron_hopping_cutoff) {
			cout << "Error! The event recalculation cutoff radius must not be less than the polaron hopping cutoff radius." << endl;
			return false;
		}
		if (Enable_selective_recalc && Recalc_cutoff < Exciton_dissociation_cutoff) {
			cout << "Error! The event recalculation cutoff radius must not be less than the exciton dissociation cutoff radius." << endl;
			return false;
		}
		// Check film architecture parameters
		if (Enable_bilayer && Thickness_donor + Thickness_acceptor != Params_lattice.Height) {
			cout << "Error! When using the bilayer film architecture, the sum of the donor and the acceptor thicknesses must equal the lattice height." << endl;
			return false;
		}
		// Possible device architectures:
		// Neat
		// Bilayer
		// Random Blend
		// Import morphology
		int N_architectures_enabled = 0;
		if (Enable_neat) {
			N_architectures_enabled++;
		}
		if (Enable_bilayer) {
			N_architectures_enabled++;
		}
		if (Enable_random_blend) {
			N_architectures_enabled++;
		}
		if (Enable_import_morphology) {
			N_architectures_enabled++;
		}
		if (N_architectures_enabled == 0) {
			cout << "Error! A film architecture must be enabled." << endl;
			return false;
		}
		if (N_architectures_enabled > 1) {
			cout << "Error! Only one film architecture can be enabled." << endl;
			return false;
		}
		// Check test parameters
		if (Enable_ToF_test && ((Enable_ToF_random_placement && Enable_ToF_energy_placement) || (!Enable_ToF_random_placement && !Enable_ToF_energy_placement))) {
			cout << "Error! For a time-of-flight charge transport test either the random placement or the low energy placement option must be enabled." << endl;
			return false;
		}
		if (Enable_ToF_test && Enable_bilayer) {
			cout << "Error! The bilayer film architecture cannot be used with the time-of-flight charge transport test." << endl;
			return false;
		}
		if (Enable_ToF_test && Params_lattice.Enable_periodic_z) {
			cout << "Error! The z-direction periodic boundary must be disabled in order to run the time-of-flight charge transport test." << endl;
			return false;
		}
		if (Enable_ToF_test && !ToF_polaron_type && Enable_neat) {
			cout << "Error! The time-of-flight charge transport test cannot be performed with electrons on a neat film architecture. Use holes instead." << endl;
			return false;
		}
		if (Enable_IQE_test && Params_lattice.Enable_periodic_z) {
			cout << "Error! The z-direction periodic boundary must be disabled in order to run the internal quantum efficiency test." << endl;
			return false;
		}
		if (Enable_neat && Enable_IQE_test) {
			cout << "Error! The neat film architecture cannot be used with the internal quantum efficiency test." << endl;
			return false;
		}
		if (!(N_tests > 0)) {
			cout << "Error! The number of tests must be greater than zero." << endl;
			return false;
		}
		// Possible simulation tests:
		// Exciton diffusion test
		// ToF test
		// IQE test
		// Dynamics test
		int N_tests_enabled = 0;
		if (Enable_exciton_diffusion_test)
			N_tests_enabled++;
		if (Enable_ToF_test)
			N_tests_enabled++;
		if (Enable_IQE_test)
			N_tests_enabled++;
		if (Enable_dynamics_test)
			N_tests_enabled++;
		if (Enable_steady_transport_test)
			N_tests_enabled++;
		if (N_tests_enabled > 1) {
			cout << "Error! Only one test can be enabled." << endl;
			return false;
		}
		if (N_tests_enabled == 0) {
			cout << "Error! One of the tests must be enabled." << endl;
			return false;
		}
		// Check dynamics test conditions
		if (Enable_dynamics_test && !Enable_dynamics_extraction && Internal_potential != 0) {
			cout << "Error! When running a dynamics test without extraction, the internal potential must be set to zero." << endl;
			return false;
		}
		if (Enable_dynamics_test && Enable_dynamics_extraction && Params_lattice.Enable_periodic_z) {
			cout << "Error! When running a dynamics test with extraction, z-direction periodic boundaries cannot be used." << endl;
			return false;
		}
		// Check steady transport test parameters
		if (Enable_steady_transport_test && !Params_lattice.Enable_periodic_z) {
			cout << "Error! When running a steady transport test, z-direction periodic boundaries must be used." << endl;
			return false;
		}
		if (Enable_steady_transport_test && abs(Internal_potential) < 1e-9) {
			cout << "Error! When running a steady transport test, the internal potential must not be set to zero." << endl;
			return false;
		}
		if (Enable_steady_transport_test && N_equilibration_events < 0) {
			cout << "Error! When running a steady transport test, the number of equilibration events cannot be negative." << endl;
			return false;
		}
		if (Enable_steady_transport_test && round_int(Steady_carrier_density*Params_lattice.Length*Params_lattice.Width*Params_lattice.Height*intpow(Params_lattice.Unit_size*1e-7,3)) < 1) {
			cout << "Error! When running a steady transport test, the steady carrier density must be large enough that there is at least one polaron in the lattice." << endl;
			return false;
		}
		if (Enable_steady_transport_test && Enable_bilayer) {
			cout << "Error! When running a steady transport test, the bilayer device architecture cannot be used." << endl;
			return false;
		}
		// Check exciton parameters
		if (Exciton_generation_rate_donor < 0 || Exciton_generation_rate_acceptor < 0) {
			cout << "Error! The exciton generation rate of the donor and acceptor must not be negative." << endl;
			return false;
		}
		if (!(Singlet_lifetime_donor > 0) || !(Singlet_lifetime_acceptor > 0)) {
			cout << "Error! The singlet exciton lifetime of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(Triplet_lifetime_donor > 0) || !(Triplet_lifetime_acceptor > 0)) {
			cout << "Error! The triplet exciton lifetime of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(R_singlet_hopping_donor > 0) || !(R_singlet_hopping_acceptor > 0)) {
			cout << "Error! The singlet exciton hopping rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(Singlet_localization_donor > 0) || !(Singlet_localization_acceptor > 0)) {
			cout << "Error! The singlet exciton localization parameter of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(R_triplet_hopping_donor > 0) || !(R_triplet_hopping_acceptor > 0)) {
			cout << "Error! The triplet exciton hopping rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(Triplet_localization_donor > 0) || !(Triplet_localization_acceptor > 0)) {
			cout << "Error! The triplet exciton localization parameter of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(R_exciton_exciton_annihilation_donor > 0) || !(R_exciton_exciton_annihilation_acceptor > 0)) {
			cout << "Error! The exciton-exciton annihilation rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(R_exciton_polaron_annihilation_donor > 0) || !(R_exciton_polaron_annihilation_acceptor > 0)) {
			cout << "Error! The exciton-polaron annihilation rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(FRET_cutoff > 0)) {
			cout << "Error! The FRET cutoff radius must be greater than zero." << endl;
			return false;
		}
		if (E_exciton_binding_donor < 0 || E_exciton_binding_acceptor < 0) {
			cout << "Error! The exciton binding energy of the donor and acceptor cannot be negative." << endl;
			return false;
		}
		if (!(R_exciton_dissociation_donor > 0) || !(R_exciton_dissociation_acceptor > 0)) {
			cout << "Error! The exciton dissociation rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(Exciton_dissociation_cutoff > 0)) {
			cout << "Error! The exciton dissociation cutoff radius must be greater than zero." << endl;
			return false;
		}
		if (!(R_exciton_isc_donor > 0) || !(R_exciton_isc_acceptor > 0)) {
			cout << "Error! The exciton intersystem crossing rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(R_exciton_risc_donor > 0) || !(R_exciton_risc_acceptor > 0)) {
			cout << "Error! The exciton reverse intersystem crossing rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (E_exciton_ST_donor < 0 || E_exciton_ST_acceptor < 0) {
			cout << "Error! The exciton singlet-triplet splitting energy of the donor and acceptor must not be negative." << endl;
			return false;
		}
		// Check polaron parameters
		if (!(R_polaron_hopping_donor > 0) || !(R_polaron_hopping_acceptor > 0)) {
			cout << "Error! The polaron hopping rate of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (!(Polaron_localization_donor > 0) || !(Polaron_localization_acceptor > 0)) {
			cout << "Error! The polaron localization parameter of the donor and acceptor must be greater than zero." << endl;
			return false;
		}
		if (Enable_miller_abrahams && Enable_marcus) {
			cout << "Error! The Miller-Abrahams and the Marcus polaron hopping models cannot both be enabled." << endl;
			return false;
		}
		if (!Enable_miller_abrahams && !Enable_marcus) {
			cout << "Error! Either the Miller-Abrahams or the Marcus polaron hopping model must be enabled." << endl;
			return false;
		}
		if (Reorganization_donor < 0 || Reorganization_acceptor < 0) {
			cout << "Error! The polaron reorganization energy of the donor and acceptor must not be negative." << endl;
			return false;
		}
		if (!(R_polaron_recombination > 0)) {
			cout << "Error! The polaron recombination rate must be greater than zero." << endl;
			return false;
		}
		if (!(Polaron_hopping_cutoff > 0)) {
			cout << "Error! The polaron hopping cutoff radius must be greater than zero." << endl;
			return false;
		}
		if (!(Polaron_delocalization_length > 0)) {
			cout << "Error! The polaron delocalization length must be greater than zero." << endl;
			return false;
		}
		// Check lattice site parameters
		if (Homo_donor < 0 || Lumo_donor < 0) {
			cout << "Error! The HOMO and LUMO parameters of the donor must not be negative." << endl;
			return false;
		}
		if (Homo_acceptor < 0 || Lumo_acceptor < 0) {
			cout << "Error! The HOMO and LUMO parameters of the acceptor must not be negative." << endl;
			return false;
		}
		if (Enable_gaussian_dos && Enable_exponential_dos) {
			cout << "Error! The Gaussian and exponential disorder models cannot both be enabled." << endl;
			return false;
		}
		if (Enable_gaussian_dos && (Energy_stdev_donor < 0 || Energy_stdev_acceptor < 0)) {
			cout << "Error! When using the Gaussian disorder model, the standard deviation cannot be negative." << endl;
			return false;
		}
		if (Enable_exponential_dos && (Energy_urbach_donor < 0 || Energy_urbach_acceptor < 0)) {
			cout << "Error! When using the exponential disorder model, the Urbach energy cannot be negative." << endl;
			return false;
		}
		int kernel_counter = 0;
		if (Enable_correlated_disorder && Enable_gaussian_kernel) {
			kernel_counter++;
		}
		if (Enable_correlated_disorder && Enable_power_kernel) {
			kernel_counter++;
		}
		if (Enable_correlated_disorder && kernel_counter != 1) {
			cout << "Error! When using the correlated disorder model, you must enable one and only one kernel. You have " << kernel_counter << " kernels enabled." << endl;
			return false;
		}
		if (Enable_correlated_disorder && Enable_power_kernel && !(Power_kernel_exponent == -1 || Power_kernel_exponent == -2)) {
			cout << "Error! When using the correlated disorder model with the power kernel, the power kernel exponent must be either -1 or -2." << endl;
			return false;
		}
		if (Enable_correlated_disorder && (Disorder_correlation_length < 0.999 || Disorder_correlation_length > 2.001)) {
			cout << "Error! When using the correlated disorder model, the disorder correlation length must be in the range between 1.0 and 2.0." << endl;
			return false;
		}
		if (Enable_correlated_disorder && !Enable_neat && Energy_stdev_donor != Energy_stdev_acceptor) {
			cout << "Error! When using the correlated disorder model, the standard deviation of the donor and acceptor Gaussian DOS must be equal." << endl;
			return false;
		}
		if (Enable_correlated_disorder && !Enable_gaussian_dos) {
			cout << "Error! The correlated disorder model can only be used with a Gaussian density of states." << endl;
			return false;
		}
		if (Enable_interfacial_energy_shift && (Energy_shift_donor < 0 || Energy_shift_acceptor < 0)) {
			cout << "Error! When enabling the interfacial energy shift model, the energy shift factor for the donor and the acceptor phases must not be negative." << endl;
			return false;
		}
		if (Enable_import_energies && (Enable_gaussian_dos || Enable_exponential_dos)) {
			cout << "Error! When importing site energies from a file, the Gaussian and exponential density of states models must not be enabled." << endl;
			return false;
		}
		if (Enable_import_energies && Enable_interfacial_energy_shift) {
			cout << "Error! When importing site energies from a file, the interfacial energy shift model must not be enabled." << endl;
			return false;
		}
		if (Enable_import_energies && (int)Energies_import_filename.size() == 0) {
			cout << "Error! When importing site energies from a file, a valid filename must be provided." << endl;
			return false;
		}
        if (Output_interval < 1) {
            cout << "Error! Output interval is smaller than 1." << endl;
            return false;
        }
        
		// Check Coulomb interaction parameters
		if (!(Coulomb_cutoff > 0)) {
			cout << "Error! The Coulomb cutoff radius must be greater than zero." << endl;
			return false;
		}
		if (!(Dielectric_donor > 0) || !(Dielectric_acceptor > 0)) {
			cout << "Error! The dielectric constant of the donor and the acceptor must be greater than zero." << endl;
			return false;
		}
		return true;
	}

	bool Parameters::importParameters(ifstream& inputfile) {
		string line;
		string var;
		size_t pos;
		vector<string> stringvars;
		bool Error_found = false;
		if (!inputfile.is_open() || !inputfile) {
			throw invalid_argument("Error importing parameter file because ifstream cannot read the parameter file.");
		}
		while (inputfile.good()) {
			getline(inputfile, line);
			if ((line.substr(0, 2)).compare("--") != 0 && (line.substr(0, 2)).compare("##") != 0 && line.compare("") != 0) {
				pos = line.find_first_of("/", 0);
				var = line.substr(0, pos);
				var = removeWhitespace(var);
				stringvars.push_back(var);
			}
		}
		int i = 0;
		// KMC Algorithm Parameters
		try {
			Enable_FRM = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting first reaction method option." << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_selective_recalc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting selective recalculation method option." << endl;
			Error_found = true;
		}
		i++;
		Recalc_cutoff = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_full_recalc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting full recalculation method option." << endl;
			Error_found = true;
		}
		i++;
		//enable_periodic_x
		try {
			Params_lattice.Enable_periodic_x = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting x-periodic boundary options" << endl;
			Error_found = true;
		}
		i++;
		//enable_periodic_y
		try {
			Params_lattice.Enable_periodic_y = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting y-periodic boundary options" << endl;
			Error_found = true;
		}
		i++;
		//enable_periodic_z
		try {
			Params_lattice.Enable_periodic_z = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting z-periodic boundary options" << endl;
			Error_found = true;
		}
		i++;
		Params_lattice.Length = atoi(stringvars[i].c_str());
		i++;
		Params_lattice.Width = atoi(stringvars[i].c_str());
		i++;
		Params_lattice.Height = atoi(stringvars[i].c_str());
		i++;
		Params_lattice.Unit_size = atof(stringvars[i].c_str());
		i++;
		Temperature = atoi(stringvars[i].c_str());
		i++;
		Internal_potential = atof(stringvars[i].c_str());
		i++;
		// Film Architecture Parameters
		try {
			Enable_neat = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling neat film architecture." << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_bilayer = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling bilayer film architecture." << endl;
			Error_found = true;
		}
		i++;
		Thickness_donor = atoi(stringvars[i].c_str());
		i++;
		Thickness_acceptor = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_random_blend = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling random blend film architecture." << endl;
			Error_found = true;
		}
		i++;
		Acceptor_conc = atof(stringvars[i].c_str());;
		i++;
		try {
			Enable_import_morphology_single = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling morphology import." << endl;
			Error_found = true;
		}
		i++;
		Morphology_filename = stringvars[i];
		i++;
		try {
			Enable_import_morphology_set = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling morphology set import." << endl;
			Error_found = true;
		}
		i++;
		Morphology_set_format = stringvars[i];
		i++;
		N_test_morphologies = atoi(stringvars[i].c_str());
		i++;
		N_morphology_set_size = atoi(stringvars[i].c_str());
		i++;
		if (Enable_import_morphology_single || Enable_import_morphology_set) {
			Enable_import_morphology = true;
		}
		else {
			Enable_import_morphology = false;
		}
		// Test Parameters
		N_tests = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_exciton_diffusion_test = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling the exciton diffusion test." << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_ToF_test = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling the time-of-flight polaron transport test." << endl;
			Error_found = true;
		}
		i++;
		if (stringvars[i].compare("electron") == 0) {
			ToF_polaron_type = false;
		}
		else if (stringvars[i].compare("hole") == 0) {
			ToF_polaron_type = true;
		}
		else {
			cout << "Error setting polaron type for the time-of-flight test." << endl;
			return false;
		}
		i++;
		ToF_initial_polarons = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_ToF_random_placement = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling ToF random placement option." << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_ToF_energy_placement = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling ToF low energy placement option." << endl;
			Error_found = true;
		}
		i++;
		ToF_placement_energy = atof(stringvars[i].c_str());
		i++;
		ToF_transient_start = atof(stringvars[i].c_str());
		i++;
		ToF_transient_end = atof(stringvars[i].c_str());
		i++;
		ToF_pnts_per_decade = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_IQE_test = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling the internal quantum efficiency test." << endl;
			Error_found = true;
		}
		i++;
		IQE_time_cutoff = atof(stringvars[i].c_str());
		i++;
		try {
			Enable_extraction_map_output = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting charge extraction map output settings." << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_dynamics_test = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling the dynamics test." << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_dynamics_extraction = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting dynamics test extraction option." << endl;
			Error_found = true;
		}
		i++;
		Dynamics_initial_exciton_conc = atof(stringvars[i].c_str());
		i++;
		Dynamics_transient_start = atof(stringvars[i].c_str());
		i++;
		Dynamics_transient_end = atof(stringvars[i].c_str());
		i++;
		Dynamics_pnts_per_decade = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_steady_transport_test = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error enabling the steady transport test." << endl;
			Error_found = true;
		}
		i++;
		Steady_carrier_density = atof(stringvars[i].c_str());
		i++;
		N_equilibration_events = atoi(stringvars[i].c_str());
		i++;
		// Exciton Parameters
		Exciton_generation_rate_donor = atof(stringvars[i].c_str());
		i++;
		Exciton_generation_rate_acceptor = atof(stringvars[i].c_str());
		i++;
		Singlet_lifetime_donor = atof(stringvars[i].c_str());
		i++;
		Singlet_lifetime_acceptor = atof(stringvars[i].c_str());
		i++;
		Triplet_lifetime_donor = atof(stringvars[i].c_str());
		i++;
		Triplet_lifetime_acceptor = atof(stringvars[i].c_str());
		i++;
		R_singlet_hopping_donor = atof(stringvars[i].c_str());
		i++;
		R_singlet_hopping_acceptor = atof(stringvars[i].c_str());
		i++;
		Singlet_localization_donor = atof(stringvars[i].c_str());
		i++;
		Singlet_localization_acceptor = atof(stringvars[i].c_str());
		i++;
		R_triplet_hopping_donor = atof(stringvars[i].c_str());
		i++;
		R_triplet_hopping_acceptor = atof(stringvars[i].c_str());
		i++;
		Triplet_localization_donor = atof(stringvars[i].c_str());
		i++;
		Triplet_localization_acceptor = atof(stringvars[i].c_str());
		i++;
		try {
			Enable_FRET_triplet_annihilation = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting FRET triplet annihilation option." << endl;
			return false;
		}
		i++;
		R_exciton_exciton_annihilation_donor = atof(stringvars[i].c_str());
		i++;
		R_exciton_exciton_annihilation_acceptor = atof(stringvars[i].c_str());
		i++;
		R_exciton_polaron_annihilation_donor = atof(stringvars[i].c_str());
		i++;
		R_exciton_polaron_annihilation_acceptor = atof(stringvars[i].c_str());
		i++;
		FRET_cutoff = atoi(stringvars[i].c_str());
		i++;
		E_exciton_binding_donor = atof(stringvars[i].c_str());
		i++;
		E_exciton_binding_acceptor = atof(stringvars[i].c_str());
		i++;
		R_exciton_dissociation_donor = atof(stringvars[i].c_str());
		i++;
		R_exciton_dissociation_acceptor = atof(stringvars[i].c_str());
		i++;
		Exciton_dissociation_cutoff = atoi(stringvars[i].c_str());
		i++;
		R_exciton_isc_donor = atof(stringvars[i].c_str());
		i++;
		R_exciton_isc_acceptor = atof(stringvars[i].c_str());
		i++;
		R_exciton_risc_donor = atof(stringvars[i].c_str());
		i++;
		R_exciton_risc_acceptor = atof(stringvars[i].c_str());
		i++;
		E_exciton_ST_donor = atof(stringvars[i].c_str());
		i++;
		E_exciton_ST_acceptor = atof(stringvars[i].c_str());
		i++;
		// Polaron Parameters
		try {
			Enable_phase_restriction = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting polaron phase restriction option." << endl;
			Error_found = true;
		}
		i++;
		R_polaron_hopping_donor = atof(stringvars[i].c_str());
		i++;
		R_polaron_hopping_acceptor = atof(stringvars[i].c_str());
		i++;
		Polaron_localization_donor = atof(stringvars[i].c_str());
		i++;
		Polaron_localization_acceptor = atof(stringvars[i].c_str());
		i++;
		try {
			Enable_miller_abrahams = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Miller-Abrahams polaron hopping model options" << endl;
			Error_found = true;
		}
		i++;
		try {
			Enable_marcus = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Marcus polaron hopping model options" << endl;
			Error_found = true;
		}
		i++;
		Reorganization_donor = atof(stringvars[i].c_str());
		i++;
		Reorganization_acceptor = atof(stringvars[i].c_str());
		i++;
		R_polaron_recombination = atof(stringvars[i].c_str());
		i++;
		Polaron_hopping_cutoff = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_gaussian_polaron_delocalization = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Gaussian polaron delocalization option." << endl;
			Error_found = true;
		}
		i++;
		Polaron_delocalization_length = atof(stringvars[i].c_str());
		i++;
		// Lattice Parameters
		Homo_donor = atof(stringvars[i].c_str());
		i++;
		Lumo_donor = atof(stringvars[i].c_str());
		i++;
		Homo_acceptor = atof(stringvars[i].c_str());
		i++;
		Lumo_acceptor = atof(stringvars[i].c_str());
		i++;
		//enable_gaussian_dos
		try {
			Enable_gaussian_dos = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Gaussian DOS options" << endl;
			Error_found = true;
		}
		i++;
		Energy_stdev_donor = atof(stringvars[i].c_str());
		i++;
		Energy_stdev_acceptor = atof(stringvars[i].c_str());
		i++;
		//enable_exponential_dos
		try {
			Enable_exponential_dos = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Exponential DOS options" << endl;
			Error_found = true;
		}
		i++;
		Energy_urbach_donor = atof(stringvars[i].c_str());
		i++;
		Energy_urbach_acceptor = atof(stringvars[i].c_str());
		i++;
		//enable_correlated_disorder
		try {
			Enable_correlated_disorder = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Correlated Disorder options" << endl;
			Error_found = true;
		}
		i++;
		Disorder_correlation_length = atof(stringvars[i].c_str());
		i++;
		//enable_gaussian_kernel
		try {
			Enable_gaussian_kernel = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Correlated Disorder Gaussian kernel options" << endl;
			Error_found = true;
		}
		i++;
		//enable_power_kernel
		try {
			Enable_power_kernel = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Correlated Disorder Gaussian kernel options" << endl;
			Error_found = true;
		}
		i++;
		Power_kernel_exponent = atoi(stringvars[i].c_str());
		i++;
		//enable_interfacial_energy_shift
		try {
			Enable_interfacial_energy_shift = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting the interfacial energy shift options" << endl;
			Error_found = true;
		}
		i++;
		Energy_shift_donor = atof(stringvars[i].c_str());
		i++;
		Energy_shift_acceptor = atof(stringvars[i].c_str());
		i++;
		//Enable_import_energies
		try {
			Enable_import_energies = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting the import site energies options" << endl;
			Error_found = true;
		}
		i++;
		Energies_import_filename = stringvars[i];
		i++;

		//Enable_export_energies
		try {
			Enable_export_energies = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting the export site energies options" << endl;
			Error_found = true;
		}
		i++;
		Energies_export_filename = stringvars[i];
		i++;

		//Enable_import_occupancies
		try {
			Enable_import_occupancies = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting the import site occupancies options" << endl;
			Error_found = true;
		}
		i++;
		Occupancies_import_filename = stringvars[i];
		i++;        
		//Enable_export_occupancies
		try {
			Enable_export_occupancies = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting the export site occupancies options" << endl;
			Error_found = true;
		}
		i++;
		try {
			Keep_only_newest_occupancy = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting the keep only newest option" << endl;
			Error_found = true;
        }
        i++;
        Output_interval = atoi(stringvars[i].c_str());
        i++;
		Occupancies_export_filename = stringvars[i];
		i++;
		// Coulomb Calculation Parameters
		Dielectric_donor = atof(stringvars[i].c_str());
		//i++;
		Dielectric_acceptor = atof(stringvars[i].c_str());
		i++;
		//enable_coulomb_maximum
		bool Enable_Coulomb_maximum = false;
		try {
			Enable_Coulomb_maximum = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Coulomb interaction options" << endl;
			Error_found = true;
		}
		i++;
		//enable_coulomb_cutoff
		bool Enable_Coulomb_cutoff = false;
		try {
			Enable_Coulomb_cutoff = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Coulomb interaction options" << endl;
			Error_found = true;
		}
		i++;
		Coulomb_cutoff = atoi(stringvars[i].c_str());
		i++;
		if (Enable_Coulomb_maximum && Enable_Coulomb_cutoff) {
			cout << "Error! Cannot enable both the maximum Coulomb cutoff and enable use of a specific cutoff distance." << endl;
			return false;
		}
		if (!Enable_Coulomb_maximum && !Enable_Coulomb_cutoff) {
			cout << "Error! One of the Coulomb interactions calculation options must be enabled." << endl;
			return false;
		}
		if (Enable_Coulomb_maximum) {
			auto vec = { Params_lattice.Length, Params_lattice.Width, Params_lattice.Height };
			Coulomb_cutoff = (int)floor(*min_element(vec.begin(), vec.end()) / 2.0);
		}
		if (Error_found) {
			return false;
		}
		return true;
	}

}
