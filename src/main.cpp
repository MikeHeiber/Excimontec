// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "OSC_Sim.h"
#include "Parameters.h"
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <functional>

using namespace std;
using namespace Excimontec;
using namespace KMC_Lattice;

int main(int argc, char *argv[]) {
	string version = "v1.0.0-rc.4";
	// Parameters
	bool End_sim = false;
	// File declaration
	ifstream parameterfile;
	ofstream logfile;
	ofstream resultsfile;
	ofstream analysisfile;
	// Initialize variables
	string logfilename;
	Parameters params;
	int nproc = 1;
	int procid = 0;
	int elapsedtime;
	time_t time_start, time_end;
	bool success;
	bool all_finished = false;
	vector<bool> proc_finished;
	vector<bool> error_status_vec;
	vector<string> error_messages;
	char error_found = (char)0;
	// Start timer
	time_start = time(NULL);
	// Check command line arguments
	if (argc < 2) {
		cout << "Error! You must input the parameter file name as a command line argument." << endl;
		return 0;
	}
	// Check for command line enabled logging
	// Set default
	params.Enable_logging = false;
	if (argc == 3) {
		string argument(argv[2]);
		if (argument.compare("-enable_logging") == 0) {
			params.Enable_logging = true;
		}
		else {
			cout << "Error! Invalid command line argument." << endl;
			return 0;
		}
	}
	// Check for too many command line arguments
	if (argc > 3) {
		cout << "Error! Too many command line arguments." << endl;
		return 0;
	}
	// Import parameters and options from parameter file and command line arguments
	cout << "Loading input parameters from file... " << endl;
	parameterfile.open(argv[1], ifstream::in);
	if (!parameterfile.good()) {
		cout << "Error loading parameter file.  Program will now exit." << endl;
		return 0;
	}
	success = params.importParameters(parameterfile);
	parameterfile.close();
	if (!success) {
		cout << "Error importing parameters from parameter file.  Program will now exit." << endl;
		return 0;
	}
	cout << "Parameter loading complete!" << endl;
	// Initialize mpi options
	cout << "Initializing MPI options... ";
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	cout << procid << ": MPI initialization complete!" << endl;
	// Initialize error monitoring vectors
	proc_finished.assign(nproc, false);
	error_status_vec.assign(nproc, false);
	error_messages.assign(nproc, "");
	// Morphology set import handling
	if (params.Enable_import_morphology_set && params.N_test_morphologies > nproc) {
		cout << "Error! The number of requested processors cannot be less than the number of morphologies tested." << endl;
		cout << "You have requested " << nproc << " processors for " << params.N_test_morphologies << " morphologies." << endl;
		return 0;
	}
	if (params.Enable_import_morphology_set && params.N_test_morphologies > params.N_morphology_set_size) {
		cout << "Error! The number of tested morphologies cannot be greater than the number of morphologies in the set." << endl;
		cout << "You have asked to test " << params.N_test_morphologies << " morphologies out of a " << params.N_morphology_set_size << " morphology set." << endl;
		return 0;
	}
	if (params.Enable_import_morphology_set) {
		int* selected_morphologies = new int[nproc];
		if (procid == 0) {
			default_random_engine generator((int)time(0));
			vector<int> morphology_set_original(params.N_test_morphologies);
			vector<int> morphology_set;
			// Select morphologies from the morphology set
			for (int n = 0; n < params.N_morphology_set_size; n++) {
				morphology_set.push_back(n);
			}
			shuffle(morphology_set.begin(), morphology_set.end(), generator);
			for (int n = 0; n < params.N_test_morphologies; n++) {
				morphology_set_original[n] = morphology_set[n];
			}
			// Assign morphologies to each processor
			morphology_set.clear();
			for (int n = 0; n < nproc; n++) {
				// Fill morphology set when empty and shuffle
				if ((int)morphology_set.size() == 0) {
					morphology_set = morphology_set_original;
					shuffle(morphology_set.begin(), morphology_set.end(), generator);
				}
				// Assign morphology from the back of the set and remove it from the set
				selected_morphologies[n] = morphology_set.back();
				morphology_set.pop_back();
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(selected_morphologies, nproc, MPI_INT, 0, MPI_COMM_WORLD);
		// Parse input morphology set file format
		int pos = (int)params.Morphology_set_format.find("#");
		string prefix = params.Morphology_set_format.substr(0, pos);
		string suffix = params.Morphology_set_format.substr(pos + 1);
		cout << procid << ": Morphology " << selected_morphologies[procid] << " selected." << endl;
		params.Morphology_filename = prefix + to_string(selected_morphologies[procid]) + suffix;
		cout << procid << ": " << params.Morphology_filename << " selected." << endl;
		// Cleanup
		delete[] selected_morphologies;
	}
	// Setup file output
	cout << procid << ": Creating output files..." << endl;
	if (params.Enable_logging) {
		logfilename = "log" + to_string(procid) + ".txt";
		logfile.open(logfilename);
	}
	params.Logfile = &logfile;
	// Initialize Simulation
	cout << procid << ": Initializing simulation " << procid << "..." << endl;
	OSC_Sim sim;
	success = sim.init(params, procid);
	if (!success) {
		cout << procid << ": Initialization failed, simulation will now terminate." << endl;
		return 0;
	}
	cout << procid << ": Simulation initialization complete" << endl;
	if (params.Enable_exciton_diffusion_test) {
		cout << procid << ": Starting exciton diffusion test..." << endl;
	}
	else if (params.Enable_dynamics_test) {
		cout << procid << ": Starting dynamics test..." << endl;
	}
	else if (params.Enable_ToF_test) {
		cout << procid << ": Starting time-of-flight charge transport test..." << endl;
	}
	else if (params.Enable_IQE_test) {
		cout << procid << ": Starting internal quantum efficiency test..." << endl;
	}
	else if (params.Enable_steady_transport_test) {
		cout << procid << ": Starting steady state charge transport test..." << endl;
	}
	// Begin Simulation loop
	// Simulation ends for all procs with procid >0 when End_sim is true
	// Proc 0 only ends when End_sim is true and all_finished is true
	while (!End_sim || (procid == 0 && !all_finished)) {
		if (!End_sim) {
			success = sim.executeNextEvent();
			if (!success) {
				cout << procid << ": Event execution failed, simulation will now terminate." << endl;
			}
			End_sim = sim.checkFinished();
		}
		// Check for errors
		if (!success || sim.getN_events_executed() % 500000 == 0 || End_sim) {
			// Send completion status, error status, message length, and message content to proc 0
			char finished_status;
			char error_status;
			char msg_length;
			for (int i = 1; i < nproc; i++) {
				// Send status messages to proc 0
				if (procid == i) {
					finished_status = End_sim ? (char)1 : (char)0;
					error_status = !success ? (char)1 : (char)0;
					MPI_Send(&error_status, 1, MPI_CHAR, 0, i, MPI_COMM_WORLD);
					// If the proc has an error, send the error message to proc 0
					if (!success) {
						msg_length = (char)sim.getErrorMessage().size();
						MPI_Send(&msg_length, 1, MPI_CHAR, 0, i, MPI_COMM_WORLD);
						char* error_msg = &sim.getErrorMessage()[0];
						MPI_Send(&error_msg, (int)msg_length, MPI_CHAR, 0, i, MPI_COMM_WORLD);
					}
					MPI_Send(&finished_status, 1, MPI_CHAR, 0, i, MPI_COMM_WORLD);
				}
				// Receive messages from any processors not previously finished
				if (procid == 0 && !proc_finished[i]) {
					// Receive error status message
					MPI_Recv(&error_status, 1, MPI_CHAR, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					error_status_vec[i] = (error_status == (char)1) ? true : false;
					// If the proc has an error, then receive the error message
					if (error_status_vec[i]) {
						MPI_Recv(&msg_length, 1, MPI_CHAR, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						char* error_msg = new char[(int)msg_length];
						MPI_Recv(error_msg, (int)msg_length, MPI_CHAR, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						error_messages[i] = string(error_msg);
						error_found = (char)1;
						delete[] error_msg;
					}
					MPI_Recv(&finished_status, 1, MPI_CHAR, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					proc_finished[i] = (finished_status == (char)1) ? true : false;
				}
			}
			// Check error status of proc 0
			if (!success) {
				error_found = (char)1;
			}
			// Send error status from proc 0 to all unfinished procs
			for (int i = 1; i < nproc; i++) {
				if (procid == 0) {
					MPI_Send(&error_found, 1, MPI_CHAR, i, i, MPI_COMM_WORLD);
				}
				if (procid == i && !proc_finished[i]) {
					MPI_Recv(&error_found, 1, MPI_CHAR, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
			// Update info for proc 0
			error_status_vec[0] = !success;
			proc_finished[0] = End_sim;
			error_messages[0] = sim.getErrorMessage();
			if (error_found == (char)1) {
				break;
			}
			// Update completion status
			if (procid == 0) {
				all_finished = true;
				for (int i = 0; i < nproc; i++) {
					if (!proc_finished[i]) {
						all_finished = false;
						break;
					}
				}
			}
		}
		// Output status
		if (sim.getN_events_executed() % 1000000 == 0) {
			sim.outputStatus();
		}
		// Reset logfile
		if (params.Enable_logging) {
			if (sim.getN_events_executed() % 1000 == 0) {
				logfile.close();
				logfile.open(logfilename);
			}
		}
	}
	if (params.Enable_logging) {
		logfile.close();
	}
	cout << procid << ": Simulation finished." << endl;
	time_end = time(NULL);
	elapsedtime = (int)difftime(time_end, time_start);
	// Output disorder correlation information if correlated disorder is enabled
	if (params.Enable_correlated_disorder) {
		auto dos_correlation_data = sim.getDOSCorrelationData();
		outputVectorToFile(dos_correlation_data, "DOS_correlation_data" + to_string(procid) + ".txt");
	}
	// Output result summary for each processor
	resultsfile.open("results" + to_string(procid) + ".txt");
	resultsfile << "Excimontec " << version << " Results:\n";
	resultsfile << "Calculation time elapsed is " << (double)elapsedtime / 60 << " minutes.\n";
	resultsfile << sim.getTime() << " seconds have been simulated.\n";
	resultsfile << sim.getN_events_executed() << " events have been executed.\n\n";
	if (!success) {
		resultsfile << "An error occurred during the simulation:" << endl;
		resultsfile << sim.getErrorMessage() << endl;
	}
	else {
		if (params.Enable_exciton_diffusion_test) {
			resultsfile << "Exciton diffusion test results:\n";
			resultsfile << sim.getN_excitons_created() << " excitons have been created.\n";
			resultsfile << "Exciton diffusion length is " << vector_avg(sim.getExcitonDiffusionData()) << " ± " << vector_stdev(sim.getExcitonDiffusionData()) << " nm.\n";
			resultsfile << "Exciton hop distance is " << vector_avg(sim.getExcitonHopLengthData()) << " ± " << vector_stdev(sim.getExcitonHopLengthData()) << " nm.\n";
			resultsfile << "Exciton lifetime is " << vector_avg(sim.getExcitonLifetimeData()) << " ± " << vector_stdev(sim.getExcitonLifetimeData()) << " s.\n";
		}
		else if (params.Enable_ToF_test) {
			resultsfile << "Time-of-flight charge transport test results:\n";
			if (!params.ToF_polaron_type) {
				resultsfile << sim.getN_electrons_collected() << " of " << sim.getN_electrons_created() << " electrons have been collected.\n";
			}
			else {
				resultsfile << sim.getN_holes_collected() << " of " << sim.getN_holes_created() << " holes have been collected.\n";
			}
			resultsfile << "Transit time is " << vector_avg(sim.getTransitTimeData()) << " ± " << vector_stdev(sim.getTransitTimeData()) << " s.\n";
			resultsfile << "Charge carrier mobility is " << vector_avg(sim.calculateMobilityData(sim.getTransitTimeData())) << " ± " << vector_stdev(sim.calculateMobilityData(sim.getTransitTimeData())) << " cm^2 V^-1 s^-1.\n";
		}
		if (params.Enable_dynamics_test) {
			resultsfile << "Dynamics test results:\n";
			resultsfile << sim.getN_excitons_created() << " initial excitons were created.\n";
		}
		if (params.Enable_IQE_test) {
			resultsfile << "Internal quantum efficiency test results:\n";
			resultsfile << sim.getN_excitons_created() << " excitons have been created.\n";
		}
		if (params.Enable_IQE_test || params.Enable_dynamics_test) {
			resultsfile << sim.getN_excitons_created((short)1) << " excitons were created on donor sites.\n";
			resultsfile << sim.getN_excitons_created((short)2) << " excitons were created on acceptor sites.\n";
			resultsfile << 100 * (double)(sim.getN_singlet_excitons_dissociated() + sim.getN_triplet_excitons_dissociated()) / (double)sim.getN_excitons_created() << "% of excitons have dissociated.\n";
			resultsfile << 100 * (double)sim.getN_singlet_excitons_recombined() / (double)sim.getN_excitons_created() << "% of excitons relaxed to the ground state as singlets.\n";
			resultsfile << 100 * (double)sim.getN_triplet_excitons_recombined() / (double)sim.getN_excitons_created() << "% of excitons relaxed to the ground state as triplets.\n";
			resultsfile << 100 * (double)sim.getN_singlet_singlet_annihilations() / (double)sim.getN_excitons_created() << "% of excitons were lost to singlet-singlet annihilation.\n";
			resultsfile << 100 * (double)sim.getN_singlet_triplet_annihilations() / (double)sim.getN_excitons_created() << "% of excitons were lost to singlet-triplet annihilation.\n";
			resultsfile << 100 * (double)sim.getN_triplet_triplet_annihilations() / (double)sim.getN_excitons_created() << "% of excitons were lost to triplet-triplet annihilation.\n";
			resultsfile << 100 * (double)sim.getN_singlet_polaron_annihilations() / (double)sim.getN_excitons_created() << "% of excitons were lost to singlet-polaron annihilation.\n";
			resultsfile << 100 * (double)sim.getN_triplet_polaron_annihilations() / (double)sim.getN_excitons_created() << "% of excitons were lost to triplet-polaron annihilation.\n";
			resultsfile << 100 * (double)sim.getN_geminate_recombinations() / (double)(sim.getN_singlet_excitons_dissociated() + sim.getN_triplet_excitons_dissociated()) << "% of photogenerated charges were lost to geminate recombination.\n";
			resultsfile << 100 * (double)sim.getN_bimolecular_recombinations() / (double)(sim.getN_singlet_excitons_dissociated() + sim.getN_triplet_excitons_dissociated()) << "% of photogenerated charges were lost to bimolecular recombination.\n";
			resultsfile << 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2 * (double)(sim.getN_singlet_excitons_dissociated() + sim.getN_triplet_excitons_dissociated())) << "% of photogenerated charges were extracted.\n";
		}
		if (params.Enable_IQE_test) {
			resultsfile << "IQE = " << 100 * (double)(sim.getN_electrons_collected() + sim.getN_holes_collected()) / (2 * (double)sim.getN_excitons_created()) << "% with an internal potential of " << params.Internal_potential << " V." << endl;
		}
		if (params.Enable_steady_transport_test) {
			resultsfile << "Steady state charge transport test results:\n";
			resultsfile << "Under the conditions:\n";
			resultsfile << "Temperature = " << sim.getTemp() << endl;
			resultsfile << "Charge carrier density = " << params.Steady_carrier_density << " cm^-3.\n";
			resultsfile << "Internal electric field = " << fabs(sim.getInternalField()) << " V cm^-1.\n\n";
			resultsfile << "Current density = " << sim.getSteadyCurrentDensity() << " mA cm^-2.\n";
			resultsfile << "Charge carrier mobility = " << sim.getSteadyMobility() << " cm^2 V^-1 s^-1.\n";
			resultsfile << "Equilibration energy (without Coulomb potential) = " << sim.getSteadyEquilibrationEnergy() << " eV.\n";
			resultsfile << "Equilibration energy (with Coulomb potential) = " << sim.getSteadyEquilibrationEnergy_Coulomb() << " eV.\n";
			resultsfile << "Transport energy (without Coulomb potential) = " << sim.getSteadyTransportEnergy() << " eV.\n";
			resultsfile << "Transport energy (with Coulomb potential) = " << sim.getSteadyTransportEnergy_Coulomb() << " eV.\n";
		}
		resultsfile << endl;
	}
	resultsfile.close();
	// Output charge extraction map data
	if (success && params.Enable_extraction_map_output && (params.Enable_ToF_test || params.Enable_IQE_test)) {
		if (params.Enable_ToF_test) {
			string filename = "Charge_extraction_map" + to_string(procid) + ".txt";
			vector<string> extraction_data = sim.getChargeExtractionMap(params.ToF_polaron_type);
			outputVectorToFile(extraction_data, filename);
		}
		if (params.Enable_IQE_test) {
			string filename = "Electron_extraction_map" + to_string(procid) + ".txt";
			vector<string> extraction_data = sim.getChargeExtractionMap(false);
			outputVectorToFile(extraction_data, filename);
			filename = "Hole_extraction_map" + to_string(procid) + ".txt";
			extraction_data = sim.getChargeExtractionMap(true);
			outputVectorToFile(extraction_data, filename);
		}
	}
	// Output overall analysis results from all processors
	int elapsedtime_sum;
	MPI_Reduce(&elapsedtime, &elapsedtime_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (procid == 0) {
		analysisfile.open("analysis_summary.txt");
		analysisfile << "Excimontec " << version << " Results Summary:\n";
		analysisfile << "Simulation was performed on " << nproc << " processors.\n";
		analysisfile << "Average calculation time was " << (double)elapsedtime_sum / (60 * nproc) << " minutes.\n\n";
		if (error_found == (char)1) {
			analysisfile << endl << "An error occurred on one or more processors:" << endl;
			for (int i = 0; i < nproc; i++) {
				if (error_status_vec[i]) {
					analysisfile << i << ": " << error_messages[i] << endl;
				}
			}
		}
	}
	if (error_found == (char)0 && params.Enable_exciton_diffusion_test) {
		vector<double> exciton_diffusion_data;
		vector<int> exciton_hop_length_data;
		vector<double> exciton_lifetime_data;
		exciton_diffusion_data = MPI_gatherVectors(sim.getExcitonDiffusionData());
		exciton_hop_length_data = MPI_gatherVectors(sim.getExcitonHopLengthData());
		exciton_lifetime_data = MPI_gatherVectors(sim.getExcitonLifetimeData());
		if (procid == 0) {
			analysisfile << "Overall exciton diffusion test results:\n";
			analysisfile << nproc * (sim.getN_singlet_excitons_recombined() + sim.getN_triplet_excitons_recombined()) << " total excitons tested." << endl;
			analysisfile << "Exciton diffusion length is " << vector_avg(exciton_diffusion_data) << " ± " << vector_stdev(exciton_diffusion_data) << " nm.\n";
			analysisfile << "Exciton hop distance is " << sqrt(vector_avg(exciton_hop_length_data))*params.Params_lattice.Unit_size << " ± " << sqrt(vector_stdev(exciton_hop_length_data))*params.Params_lattice.Unit_size << " nm.\n";
			analysisfile << "Exciton lifetime is " << vector_avg(exciton_lifetime_data) << " ± " << vector_stdev(exciton_lifetime_data) << " s.\n";
		}
	}
	if (error_found == (char)0 && params.Enable_ToF_test) {
		int N_transient_cycles = sim.getN_transient_cycles();
		int N_transient_cycles_sum;
		MPI_Reduce(&N_transient_cycles, &N_transient_cycles_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		vector<double> transit_times_all = MPI_gatherVectors(sim.getTransitTimeData());
		int transit_attempts = ((sim.getN_electrons_collected() > sim.getN_holes_collected()) ? sim.getN_electrons_created() : (sim.getN_holes_created()));
		int transit_attempts_total;
		MPI_Reduce(&transit_attempts, &transit_attempts_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		vector<int> counts = MPI_calculateVectorSum(sim.getToFTransientCounts());
		vector<double> energies = MPI_calculateVectorSum(sim.getToFTransientEnergies());
		vector<double> velocities = MPI_calculateVectorSum(sim.getToFTransientVelocities());
		vector<double> times = sim.getToFTransientTimes();
		if (procid == 0) {
			// ToF main results output
			vector<double> mobility_data_all = sim.calculateMobilityData(transit_times_all);
			double electric_field = fabs(sim.getInternalField());
			ofstream tof_resultsfile;
			tof_resultsfile.open("ToF_results.txt");
			tof_resultsfile << "Electric Field (V/cm),Transit Time Avg (s),Transit Time Stdev (s),Mobility Avg (cm^2 V^-1 s^-1),Mobility Stdev (cm^2 V^-1 s^-1)" << endl;
			tof_resultsfile << electric_field << "," << vector_avg(transit_times_all) << "," << vector_stdev(transit_times_all) << "," << vector_avg(mobility_data_all) << "," << vector_stdev(mobility_data_all) << endl;
			tof_resultsfile.close();
			// ToF transient output
			ofstream transientfile;
			transientfile.open("ToF_average_transients.txt");
			transientfile << "Time (s),Current Density (mA cm^-2),Average Mobility (cm^2 V^-1 s^-1),Average Energy (eV),Carrier Density (cm^-3)" << endl;
			double volume_total = N_transient_cycles_sum * sim.getVolume();
			for (int i = 0; i < (int)velocities.size(); i++) {
				if ((double)counts[i] > 0.95*counts[0]) {
					transientfile << times[i] << "," << 1000.0 * Elementary_charge*velocities[i] / volume_total << "," << (velocities[i] / (double)counts[i]) / electric_field << "," << energies[i] / (double)counts[i] << "," << (double)counts[i] / volume_total << endl;
				}
				else if (counts[i] > 0) {
					transientfile << times[i] << "," << 1000.0 * Elementary_charge*velocities[i] / volume_total << "," << "NaN" << "," << "NaN" << "," << (double)counts[i] / volume_total << endl;
				}
				else {
					transientfile << times[i] << ",0,NaN,NaN,0" << endl;
				}
			}
			transientfile.close();
			// ToF transit time distribution output
			ofstream transitdistfile;
			transitdistfile.open("ToF_transit_time_hist.txt");
			auto transit_dist = sim.calculateTransitTimeHist(transit_times_all, transit_attempts_total);
			transitdistfile << "Transit Time (s),Probability" << endl;
			for (int i = 0; i < (int)transit_dist.size(); i++) {
				transitdistfile << transit_dist[i].first << "," << transit_dist[i].second << endl;
			}
			transitdistfile.close();
			// Analysis Output
			analysisfile << "Overall time-of-flight charge transport test results:\n";
			if (!params.ToF_polaron_type) {
				analysisfile << nproc * sim.getN_electrons_collected() << " total electrons collected out of " << transit_attempts_total << " total attempts.\n";
			}
			else {
				analysisfile << nproc * sim.getN_holes_collected() << " total holes collected out of " << transit_attempts_total << " total attempts.\n";
			}

			analysisfile << "Transit time is " << vector_avg(transit_times_all) << " ± " << vector_stdev(transit_times_all) << " s.\n";
			analysisfile << "Charge carrier mobility is " << vector_avg(mobility_data_all) << " ± " << vector_stdev(mobility_data_all) << " cm^2 V^-1 s^-1.\n";
		}
	}
	if (error_found == (char)0 && params.Enable_dynamics_test) {
		int N_transient_cycles = sim.getN_transient_cycles();
		int N_transient_cycles_sum;
		MPI_Reduce(&N_transient_cycles, &N_transient_cycles_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		vector<double> times = sim.getDynamicsTransientTimes();
		vector<int> singlets_total = MPI_calculateVectorSum(sim.getDynamicsTransientSinglets());
		vector<int> triplets_total = MPI_calculateVectorSum(sim.getDynamicsTransientTriplets());
		vector<int> electrons_total = MPI_calculateVectorSum(sim.getDynamicsTransientElectrons());
		vector<int> holes_total = MPI_calculateVectorSum(sim.getDynamicsTransientHoles());
		vector<double> exciton_energies = MPI_calculateVectorSum(sim.getDynamicsExcitonEnergies());
		vector<double> electron_energies = MPI_calculateVectorSum(sim.getDynamicsElectronEnergies());
		vector<double> hole_energies = MPI_calculateVectorSum(sim.getDynamicsHoleEnergies());
		vector<double> exciton_msdv = MPI_calculateVectorSum(sim.getDynamicsExcitonMSDV());
		vector<double> electron_msdv = MPI_calculateVectorSum(sim.getDynamicsElectronMSDV());
		vector<double> hole_msdv = MPI_calculateVectorSum(sim.getDynamicsHoleMSDV());
		if (procid == 0) {
			ofstream transientfile;
			transientfile.open("dynamics_average_transients.txt");
			transientfile << "Time (s),Singlet Exciton Density (cm^-3),Triplet Exciton Density (cm^-3),Electron Density (cm^-3),Hole Density (cm^-3)";
			transientfile << ",Average Exciton Energy (eV),Exciton MSDV (cm^2 s^-1)";
			transientfile << ",Average Electron Energy (eV),Electron MSDV (cm^2 s^-1)";
			transientfile << ",Average Hole Energy (eV),Hole MSDV (cm^2 s^-1)" << endl;
			double volume_total = N_transient_cycles_sum * sim.getVolume();
			for (int i = 0; i < (int)times.size(); i++) {
				transientfile << times[i] << "," << singlets_total[i] / volume_total << "," << triplets_total[i] / volume_total << "," << electrons_total[i] / volume_total << "," << holes_total[i] / volume_total;
				if ((singlets_total[i] + triplets_total[i]) > 0 && (singlets_total[i] + triplets_total[i]) > 5 * N_transient_cycles_sum) {
					transientfile << "," << exciton_energies[i] / (singlets_total[i] + triplets_total[i]) << "," << exciton_msdv[i] / (singlets_total[i] + triplets_total[i]);
				}
				else if ((singlets_total[i] + triplets_total[i]) > 0) {
					transientfile << "," << "NaN" << "," << exciton_msdv[i] / (singlets_total[i] + triplets_total[i]);
				}
				else {
					transientfile << ",NaN,NaN";
				}
				if (electrons_total[i] > 0 && electrons_total[i] > 5 * N_transient_cycles_sum) {
					transientfile << "," << electron_energies[i] / electrons_total[i] << "," << electron_msdv[i] / electrons_total[i];
				}
				else if (electrons_total[i] > 0) {
					transientfile << "," << "NaN" << "," << electron_msdv[i] / electrons_total[i];
				}
				else {
					transientfile << ",NaN,NaN";
				}
				if (holes_total[i] > 0 && holes_total[i] > 5 * N_transient_cycles_sum) {
					transientfile << "," << hole_energies[i] / holes_total[i] << "," << hole_msdv[i] / holes_total[i] << endl;
				}
				else if (holes_total[i] > 0) {
					transientfile << "," << "NaN" << "," << hole_msdv[i] / holes_total[i] << endl;
				}
				else {
					transientfile << ",NaN,NaN" << endl;
				}
			}
			transientfile.close();
		}
	}
	if (error_found == (char)0 && (params.Enable_dynamics_test || params.Enable_IQE_test || params.Enable_exciton_diffusion_test)) {
		int excitons_created = sim.getN_excitons_created();
		int excitons_created_total;
		MPI_Reduce(&excitons_created, &excitons_created_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int excitons_created_donor = sim.getN_excitons_created((short)1);
		int excitons_created_donor_total;
		MPI_Reduce(&excitons_created_donor, &excitons_created_donor_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int excitons_created_acceptor = sim.getN_excitons_created((short)2);
		int excitons_created_acceptor_total;
		MPI_Reduce(&excitons_created_acceptor, &excitons_created_acceptor_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int excitons_dissociated = sim.getN_singlet_excitons_dissociated() + sim.getN_triplet_excitons_dissociated();
		int excitons_dissociated_total;
		MPI_Reduce(&excitons_dissociated, &excitons_dissociated_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int singlet_excitons_recombined = sim.getN_singlet_excitons_recombined();
		int singlet_excitons_recombined_total;
		MPI_Reduce(&singlet_excitons_recombined, &singlet_excitons_recombined_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int triplet_excitons_recombined = sim.getN_triplet_excitons_recombined();
		int triplet_excitons_recombined_total;
		MPI_Reduce(&triplet_excitons_recombined, &triplet_excitons_recombined_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int singlet_singlet_annihilations = sim.getN_singlet_singlet_annihilations();
		int singlet_singlet_annihilations_total;
		MPI_Reduce(&singlet_singlet_annihilations, &singlet_singlet_annihilations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int singlet_triplet_annihilations = sim.getN_singlet_triplet_annihilations();
		int singlet_triplet_annihilations_total;
		MPI_Reduce(&singlet_triplet_annihilations, &singlet_triplet_annihilations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int triplet_triplet_annihilations = sim.getN_triplet_triplet_annihilations();
		int triplet_triplet_annihilations_total;
		MPI_Reduce(&triplet_triplet_annihilations, &triplet_triplet_annihilations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int singlet_polaron_annihilations = sim.getN_singlet_polaron_annihilations();
		int singlet_polaron_annihilations_total;
		MPI_Reduce(&singlet_polaron_annihilations, &singlet_polaron_annihilations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int triplet_polaron_annihilations = sim.getN_triplet_polaron_annihilations();
		int triplet_polaron_annihilations_total;
		MPI_Reduce(&triplet_polaron_annihilations, &triplet_polaron_annihilations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int geminate_recombinations = sim.getN_geminate_recombinations();
		int geminate_recombinations_total;
		MPI_Reduce(&geminate_recombinations, &geminate_recombinations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int bimolecular_recombinations = sim.getN_bimolecular_recombinations();
		int bimolecular_recombinations_total;
		MPI_Reduce(&bimolecular_recombinations, &bimolecular_recombinations_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int electrons_collected = sim.getN_electrons_collected();
		int electrons_collected_total;
		MPI_Reduce(&electrons_collected, &electrons_collected_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		int holes_collected = sim.getN_holes_collected();
		int holes_collected_total;
		MPI_Reduce(&holes_collected, &holes_collected_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (procid == 0 && params.Enable_dynamics_test) {
			analysisfile << "Overall dynamics test results:\n";
		}
		if (procid == 0 && params.Enable_IQE_test) {
			analysisfile << "Overall internal quantum efficiency test results:\n";
		}
		if (procid == 0 && params.Enable_exciton_diffusion_test) {
			analysisfile << "Overall exciton mechanism statistics:\n";
		}
		if (procid == 0) {
			analysisfile << excitons_created_total << " total excitons have been created.\n";
			analysisfile << excitons_created_donor_total << " excitons were created on donor sites.\n";
			analysisfile << excitons_created_acceptor_total << " excitons were created on acceptor sites.\n";
			analysisfile << 100 * (double)excitons_dissociated_total / (double)excitons_created_total << "% of total excitons have dissociated.\n";
			analysisfile << 100 * (double)singlet_excitons_recombined_total / (double)excitons_created_total << "% of total excitons relaxed to the ground state as singlets.\n";
			analysisfile << 100 * (double)triplet_excitons_recombined_total / (double)excitons_created_total << "% of total excitons relaxed to the ground state as triplets.\n";
			analysisfile << 100 * (double)singlet_singlet_annihilations_total / (double)excitons_created_total << "% of total excitons were lost to singlet-singlet annihilation.\n";
			analysisfile << 100 * (double)singlet_triplet_annihilations_total / (double)excitons_created_total << "% of total excitons were lost to singlet-triplet annihilation.\n";
			analysisfile << 100 * (double)triplet_triplet_annihilations_total / (double)excitons_created_total << "% of total excitons were lost to triplet-triplet annihilation.\n";
			analysisfile << 100 * (double)singlet_polaron_annihilations_total / (double)excitons_created_total << "% of total excitons were lost to singlet-polaron annihilation.\n";
			analysisfile << 100 * (double)triplet_polaron_annihilations_total / (double)excitons_created_total << "% of total excitons were lost to triplet-polaron annihilation.\n";
			if (excitons_dissociated_total > 0) {
				analysisfile << 100 * (double)geminate_recombinations_total / (double)excitons_dissociated_total << "% of total photogenerated charges were lost to geminate recombination.\n";
				analysisfile << 100 * (double)bimolecular_recombinations_total / (double)excitons_dissociated_total << "% of total photogenerated charges were lost to bimolecular recombination.\n";
				analysisfile << 100 * (double)(electrons_collected_total + holes_collected_total) / (2 * (double)excitons_dissociated_total) << "% of total photogenerated charges were extracted.\n";
			}
		}
		if (procid == 0 && params.Enable_IQE_test) {
			analysisfile << "IQE = " << 100 * (double)(electrons_collected_total + holes_collected_total) / (2 * (double)excitons_created_total) << "% with an internal potential of " << params.Internal_potential << " V." << endl;
		}
	}
	if (error_found == (char)0 && params.Enable_steady_transport_test) {
		// Calculate the average DOOS
		auto doos_avg1 = MPI_calculatePairVectorAvg(sim.getSteadyDOOS());
		auto dos_avg1 = MPI_calculatePairVectorAvg(sim.getSteadyDOS());
		auto doos_avg2 = MPI_calculatePairVectorAvg(sim.getSteadyDOOS_Coulomb());
		auto dos_avg2 = MPI_calculatePairVectorAvg(sim.getSteadyDOS_Coulomb());
		// Gather results from all procs
		auto current_densities = MPI_gatherValues(sim.getSteadyCurrentDensity());
		auto mobilities = MPI_gatherValues(sim.getSteadyMobility());
		auto equilibration_energies1 = MPI_gatherValues(sim.getSteadyEquilibrationEnergy());
		auto equilibration_energies2 = MPI_gatherValues(sim.getSteadyEquilibrationEnergy_Coulomb());
		auto transport_energies1 = MPI_gatherValues(sim.getSteadyTransportEnergy());
		auto transport_energies2 = MPI_gatherValues(sim.getSteadyTransportEnergy_Coulomb());
		// Output overall results from all procs
		if (procid == 0) {
			// Output the DOOS and DOS data
			ofstream doos_file1("DOOS_data.txt");
			doos_file1 << "Energy (eV),Density (cm^-3 eV^-1)\n";
			for (auto& item : doos_avg1) {
				doos_file1 << item.first << "," << item.second << "\n";
			}
			doos_file1.close();
			ofstream doos_file2("DOOS_Coulomb_data.txt");
			doos_file2 << "Energy (eV),Density (cm^-3 eV^-1)\n";
			for (auto& item : doos_avg2) {
				doos_file2 << item.first << "," << item.second << "\n";
			}
			doos_file2.close();
			ofstream dos_file1("DOS_data.txt");
			dos_file1 << "Energy (eV),Density (cm^-3 eV^-1)\n";
			for (auto& item : dos_avg1) {
				dos_file1 << item.first << "," << item.second << "\n";
			}
			dos_file1.close();
			ofstream dos_file2("DOS_Coulomb_data.txt");
			dos_file2 << "Energy (eV),Density (cm^-3 eV^-1)\n";
			for (auto& item : dos_avg2) {
				dos_file2 << item.first << "," << item.second << "\n";
			}
			dos_file2.close();
			// Output analysis file data
			analysisfile << "Overall steady state charge transport test results:\n";
			analysisfile << "Under the conditions:\n";
			analysisfile << "Temperature = " << sim.getTemp() << " K.\n";
			analysisfile << "Charge carrier density = " << params.Steady_carrier_density << " cm^-3.\n";
			analysisfile << "Electric field = " << fabs(sim.getInternalField()) << " V cm^-1.\n\n";
			analysisfile << "Current density = " << vector_avg(current_densities) << " ± " << vector_stdev(current_densities) << " mA cm^-2.\n";
			analysisfile << "Charge carrier mobility = " << vector_avg(mobilities) << " ± " << vector_stdev(mobilities) << " cm^2 V^-1 s^-1.\n";
			analysisfile << "Equilibration energy (without Coulomb potential) = " << vector_avg(equilibration_energies1) << " ± " << vector_stdev(equilibration_energies1) << " eV.\n";
			analysisfile << "Equilibration energy (with Coulomb potential) = " << vector_avg(equilibration_energies2) << " ± " << vector_stdev(equilibration_energies2) << " eV.\n";
			analysisfile << "Transport energy (without Coulomb potential) = " << vector_avg(transport_energies1) << " ± " << vector_stdev(transport_energies1) << " eV.\n";
			analysisfile << "Transport energy (with Coulomb potential) = " << vector_avg(transport_energies2) << " ± " << vector_stdev(transport_energies2) << " eV.\n\n";
			analysisfile << "CSV formatted results:\n";
			analysisfile << "Temperature (K),Charge Carrier Density (cm^-3),Electric Field (V cm^-1),";
			analysisfile << "Current Density Avg. (mA cm^-2),Current Density Stdev. (mA cm^-2),Mobility Avg. (cm^2 V^-1 cm^-1),Mobility Stdev. (cm^2 V^-1 cm^-1),";
			analysisfile << "Equilibration Energy Avg. w/o Coulomb,Equilibration Energy Stdev. w/o Coulomb (eV),Equilibration Energy Avg. w/ Coulomb (eV),Equilibration Energy Stdev. w/ Coulomb (eV),";
			analysisfile << "Transport Energy Avg. w/o Coulomb (eV),Transport Energy Stdev. w/o Coulomb (eV),Transport Energy Avg. w/ Coulomb (eV),Transport Energy Stdev. w/ Coulomb (eV)\n";
			analysisfile << sim.getTemp() << "," << params.Steady_carrier_density << "," << fabs(sim.getInternalField()) << ",";
			analysisfile << vector_avg(current_densities) << "," << vector_stdev(current_densities) << "," << vector_avg(mobilities) << "," << vector_stdev(mobilities) << ",";
			analysisfile << vector_avg(equilibration_energies1) << "," << vector_stdev(equilibration_energies1) << "," << vector_avg(equilibration_energies2) << "," << vector_stdev(equilibration_energies2) << ",";
			analysisfile << vector_avg(transport_energies1) << "," << vector_stdev(transport_energies1) << "," << vector_avg(transport_energies2) << "," << vector_stdev(transport_energies2) << endl;
		}
	}
	if (procid == 0) {
		analysisfile.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
