// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "OSC_Sim.h"
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <functional>

using namespace std;
using namespace Utils;

struct Parameters_main{
    bool Enable_mpi;
    bool Enable_import_morphology_single;
    string Morphology_filename;
    bool Enable_import_morphology_set;
    string Morphology_set_format;
    int N_test_morphologies;
    int N_morphology_set_size;
};

//Declare Functions
bool importParameters(ifstream * inputfile,Parameters_main& params_main,Parameters_OPV& params);

int main(int argc,char *argv[]){
    string version = "v0.4-alpha";
    // Parameters
    bool End_sim = false;
    // File declaration
    ifstream parameterfile;
    ifstream morphologyfile;
    ofstream logfile;
    ofstream resultsfile;
    ofstream analysisfile;
    stringstream ss;
    // Initialize variables
    string parameterfilename;
    string logfilename;
    Parameters_main params_main;
    Parameters_OPV params_opv;
    int nproc = 1;
    int procid = 0;
    int elapsedtime;
    time_t time_start,time_end;
    bool success;
    // Start timer
    time_start = time(NULL);
    // Import parameters and options from file and command line arguments
    cout << "Loading input parameters from file... " << endl;
    parameterfilename = argv[1];
    parameterfile.open(parameterfilename.c_str(),ifstream::in);
    if(!parameterfile.good()){
        cout << "Error loading parameter file.  Program will now exit." << endl;
        return 0;
    }
    success = importParameters(&parameterfile,params_main,params_opv);
    parameterfile.close();
    if(!success){
        cout << "Error importing parameters from parameter file.  Program will now exit." << endl;
        return 0;
    }
    cout << "Parameter loading complete!" << endl;
    // Initialize mpi options
	if (params_main.Enable_mpi) {
		cout << "Initializing MPI options... ";
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
		MPI_Comm_rank(MPI_COMM_WORLD, &procid);
		cout << procid << ": MPI initialization complete!" << endl;
	}
    // Morphology set import handling
    if(params_main.Enable_mpi && params_main.Enable_import_morphology_set && params_main.N_test_morphologies>nproc){
        cout << "Error! The number of requested processors cannot be less than the number of morphologies tested." << endl;
        cout << "You have requested " << nproc << " processors for " << params_main.N_test_morphologies << " morphologies." << endl;
        return 0;
    }
	if (params_main.Enable_import_morphology_set && params_main.Enable_mpi) {
		int* selected_morphologies = (int *)malloc(sizeof(int)*nproc);
		if (procid == 0) {
			default_random_engine gen((int)time(0));
			vector<int> morphology_set_original(params_main.N_test_morphologies);
			vector<int> morphology_set;
			// Select morphologies from the morphology set
			for (int n = 0; n < params_main.N_morphology_set_size; n++) {
				morphology_set.push_back(n);
			}
			shuffle(morphology_set.begin(), morphology_set.end(), gen);
			for (int n = 0; n < params_main.N_test_morphologies; n++) {
				morphology_set_original[n] = morphology_set[n];
			}
			// Assign morphologies to each processor
			morphology_set.clear();
			for (int n = 0; n < nproc; n++) {
				// Fill morphology set when empty and shuffle
				if ((int)morphology_set.size() == 0) {
					morphology_set = morphology_set_original;
					shuffle(morphology_set.begin(), morphology_set.end(), gen);
				}
				// Assign morphology from the back of the set and remove it from the set
				selected_morphologies[n] = morphology_set.back();
				morphology_set.pop_back();
			}
		}
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(selected_morphologies,nproc,MPI_INT,0,MPI_COMM_WORLD);
        // Parse input morphology set file format
        int pos = (int)params_main.Morphology_set_format.find("#");
        string prefix = params_main.Morphology_set_format.substr(0,pos);
        string suffix = params_main.Morphology_set_format.substr(pos+1);
        cout << procid << ": Morphology " << selected_morphologies[procid] << " selected." << endl;
        ss << prefix << selected_morphologies[procid] << suffix;
        cout << procid << ": " << ss.str() << " selected." << endl;
        params_main.Morphology_filename = ss.str();
        ss.str("");
    }
    if(params_main.Enable_import_morphology_single || params_main.Enable_import_morphology_set){
        params_opv.Enable_import_morphology = true;
        morphologyfile.open(params_main.Morphology_filename.c_str(),ifstream::in);
        if(morphologyfile.good()){
            params_opv.Morphology_file = &morphologyfile;
        }
        else{
            cout << procid << ": Error opening morphology file for importing." << endl;
            return 0;
        }
    }
	else {
		params_opv.Enable_import_morphology = false;
	}
    // Setup file output
    cout << procid << ": Creating output files..." << endl;
    if(params_opv.Enable_logging){
        ss << "log" << procid << ".txt";
        logfilename = ss.str();
        logfile.open(ss.str().c_str());
        ss.str("");
    }
    params_opv.Logfile = &logfile;
    // Initialize Simulation
    cout << procid << ": Initializing simulation " << procid << "..." << endl;
	OSC_Sim sim;
	success = sim.init(params_opv, procid);
	if (!success) {
		cout << procid << ": Initialization failed, simulation will now terminate." << endl;
		return 0;
	}
    cout << procid << ": Simulation initialization complete" << endl;
    if(params_opv.Enable_exciton_diffusion_test){
        cout << procid << ": Starting exciton diffusion test..." << endl;
    }
    else if(params_opv.Enable_dynamics_test){
        cout << procid << ": Starting dynamics test..." << endl;
    }
    else if(params_opv.Enable_ToF_test){
        cout << procid << ": Starting time-of-flight charge transport test..." << endl;
    }
    else if(params_opv.Enable_IQE_test){
        cout << procid << ": Starting internal quantum efficiency test..." << endl;
    }
    // Begin Simulation loop
    while(!End_sim){
        success = sim.executeNextEvent();
        if(!success){
            cout << procid << ": Event execution failed, simulation will now terminate." << endl;
			return 0;
        }
        // Check if simulation has finished
        End_sim = sim.checkFinished();
        // Output status
        if(sim.getN_events_executed()%1000000==0){
            sim.outputStatus();
        }
        // Reset logfile
        if(params_opv.Enable_logging){
            if(sim.getN_events_executed()%10000==0){
                logfile.close();
                logfile.open(logfilename.c_str());
            }
        }
    }
    if(params_opv.Enable_logging){
        logfile.close();
    }
	cout << procid << ": Simulation finished." << endl;
	time_end = time(NULL);
	elapsedtime = (int)difftime(time_end, time_start);
	// Output results
	if (params_opv.Enable_ToF_test || params_opv.Enable_IQE_test) {
		if (params_opv.Enable_ToF_test) {
			ss << "Charge_extraction_map" << procid << ".txt";
			string filename = ss.str();
			ss.str("");
			vector<string> extraction_data = sim.getChargeExtractionMap(params_opv.ToF_polaron_type);
			outputVectorToFile(extraction_data, filename);
		}
		if (params_opv.Enable_IQE_test) {
			ss << "Electron_extraction_map" << procid << ".txt";
			string filename = ss.str();
			ss.str("");
			vector<string> extraction_data = sim.getChargeExtractionMap(false);
			outputVectorToFile(extraction_data, filename);
			ss << "Hole_extraction_map" << procid << ".txt";
			filename = ss.str();
			ss.str("");
			extraction_data = sim.getChargeExtractionMap(true);
			outputVectorToFile(extraction_data, filename);
		}
	}
    // Output result summary for each processor
    ss << "results" << procid << ".txt";
    resultsfile.open(ss.str().c_str());
    ss.str("");
    resultsfile << "Excimontec " << version << " Results:\n";
    resultsfile << "Calculation time elapsed is " << (double)elapsedtime/60 << " minutes.\n";
    resultsfile << sim.getTime() << " seconds have been simulated.\n";
    resultsfile << sim.getN_events_executed() << " events have been executed.\n";
    if(params_opv.Enable_exciton_diffusion_test){
        resultsfile << "Exciton diffusion test results:\n";
        resultsfile << sim.getN_excitons_created() << " excitons have been created.\n";
        resultsfile << "Exciton Diffusion Length is " << sim.calculateDiffusionLength_avg() << " ± " << sim.calculateDiffusionLength_stdev() << " nm.\n";
    }
    else if(params_opv.Enable_ToF_test){
        resultsfile << "Time-of-flight charge transport test results:\n";
        if(!params_opv.ToF_polaron_type){
            resultsfile << sim.getN_electrons_collected() << " of " << sim.getN_electrons_created() << " electrons have been collected.\n";
        }
        else {
            resultsfile << sim.getN_holes_collected() << " of " << sim.getN_holes_created() << " holes have been collected.\n";
        }
        resultsfile << "Transit time is " << sim.calculateTransitTime_avg() << " ± " << sim.calculateTransitTime_stdev() << " s.\n";
        resultsfile << "Charge carrier mobility is " << sim.calculateMobility_avg() << " ± " << sim.calculateMobility_stdev() << " cm^2 V^-1 s^-1.\n";
    }
    if(params_opv.Enable_dynamics_test){
        resultsfile << "Dynamics test results:\n";
        resultsfile << sim.getN_excitons_created() << " initial excitons were created.\n";
    }
    if(params_opv.Enable_IQE_test){
        resultsfile << "Internal quantum efficiency test results:\n";
        resultsfile << sim.getN_excitons_created() << " excitons have been created.\n";
    }
    if(params_opv.Enable_IQE_test || params_opv.Enable_dynamics_test){
        resultsfile << sim.getN_excitons_created((short)1) << " excitons were created on donor sites.\n";
        resultsfile << sim.getN_excitons_created((short)2) << " excitons were created on acceptor sites.\n";
        resultsfile << 100*(double)sim.getN_excitons_dissociated()/(double)sim.getN_excitons_created() << "% of excitons have dissociated.\n";
        resultsfile << 100*(double)sim.getN_geminate_recombinations()/(double)sim.getN_excitons_dissociated() << "% of photogenerated charges were lost to geminate recombination.\n";
        resultsfile << 100*(double)sim.getN_bimolecular_recombinations()/(double)sim.getN_excitons_dissociated() << "% of photogenerated charges were lost to bimolecular recombination.\n";
        resultsfile << 100*(double)(sim.getN_electrons_collected()+sim.getN_holes_collected())/(2*(double)sim.getN_excitons_dissociated()) << "% of photogenerated charges were extracted.\n";
    }
    if(params_opv.Enable_IQE_test){
        resultsfile << "IQE = " << 100*(double)(sim.getN_electrons_collected()+sim.getN_holes_collected())/(2*(double)sim.getN_excitons_created()) << "%." << endl;
    }
    resultsfile << endl;
    resultsfile.close();
    // Output overall analysis results from all processors
    if(params_main.Enable_mpi){
        int elapsedtime_sum;
        MPI_Reduce(&elapsedtime,&elapsedtime_sum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
        if(procid==0){
            ss << "analysis_summary.txt";
            analysisfile.open(ss.str().c_str());
            ss.str("");
            analysisfile << "Excimontec " << version << " Results Summary:\n";
            analysisfile << "Simulation was performed on " << nproc << " processors.\n";
            analysisfile << "Average calculation time was " << (double)elapsedtime_sum/(60*nproc) << " minutes.\n";
        }
        if(params_opv.Enable_exciton_diffusion_test){
            vector<double> diffusion_data;
            diffusion_data = MPI_gatherVectors(sim.getDiffusionData());
            if(procid==0){
                analysisfile << "Overall exciton diffusion test results:\n";
                analysisfile << nproc*sim.getN_excitons_recombined() << " total excitons tested." << endl;
                analysisfile << "Exciton diffusion length is " << vector_avg(diffusion_data) << " ± " << vector_stdev(diffusion_data) << " nm.\n";
            }
        }
        if(params_opv.Enable_ToF_test){
			int N_transient_cycles = sim.getN_transient_cycles();
			int N_transient_cycles_sum;
			MPI_Reduce(&N_transient_cycles, &N_transient_cycles_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            vector<double> transit_times = MPI_gatherVectors(sim.getTransitTimeData());
            int transit_attempts = ((sim.getN_electrons_collected()>sim.getN_holes_collected()) ? sim.getN_electrons_created() : (sim.getN_holes_created()));
            int transit_attempts_total;
            MPI_Reduce(&transit_attempts,&transit_attempts_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            vector<int> counts = MPI_calculateVectorSum(sim.getToFTransientCounts());
            vector<double> energies = MPI_calculateVectorSum(sim.getToFTransientEnergies());
            vector<double> velocities = MPI_calculateVectorSum(sim.getToFTransientVelocities());
            vector<double> times = sim.getToFTransientTimes();
            if(procid==0){
				// ToF mian results output
                vector<double> mobilities = sim.calculateMobilities(transit_times);
				double electric_field = fabs(params_opv.Bias) / (1e-7*params_opv.Height*params_opv.Unit_size);
				ofstream tof_resultsfile;
				ss << "ToF_results.txt";
				tof_resultsfile.open(ss.str().c_str());
				ss.str("");
				tof_resultsfile << "Electric Field (V/cm),Mean Transit Time (s),Mean Mobility (cm^2 V^-1 s^-1)" << endl;
				tof_resultsfile << electric_field << "," << vector_avg(transit_times) << "," << vector_avg(mobilities) << endl;
				tof_resultsfile.close();
                // ToF transient output
                ofstream transientfile;
                ss << "ToF_average_transients.txt";
                transientfile.open(ss.str().c_str());
                ss.str("");
                transientfile << "Time (s),Current (mA cm^-2),Average Mobility (cm^2 V^-1 s^-1),Average Energy (eV),Carrier Density (cm^-3)" << endl;
                double total_volume = nproc*params_opv.Length*params_opv.Width*params_opv.Height*intpow(1e-7*params_opv.Unit_size,3);
                for(int i=0;i<(int)velocities.size();i++){
                    if(counts[i]!=0){
                        transientfile << times[i] << "," << 1000*Elementary_charge*1e-7*velocities[i]/(N_transient_cycles_sum*total_volume) << "," << (velocities[i]/counts[i])/electric_field << "," << energies[i]/counts[i] << "," << counts[i]/(N_transient_cycles_sum*total_volume) << endl;
                    }
                    else{
                        transientfile << times[i] << "," << 0 << "," << 0 << "," << 0 << endl;
                    }
                }
                transientfile.close();
                // ToF transit time distribution output
                ofstream transitdistfile;
                ss << "ToF_transit_time_dist.txt";
                transitdistfile.open(ss.str().c_str());
                ss.str("");
                vector<double> transit_dist = sim.calculateTransitTimeDist(transit_times,transit_attempts_total);
                transitdistfile << "Transit Time (s),Probability" << endl;
                for(int i=0;i<(int)transit_dist.size();i++){
                    transitdistfile << times[i] << "," << transit_dist[i] << endl;
                }
                transitdistfile.close();
                // Analysis Output
                if(!params_opv.ToF_polaron_type){
                    analysisfile << nproc*sim.getN_electrons_collected() << " total electrons collected out of " << transit_attempts_total << " total attempts.\n";
                }
                else{
                    analysisfile << nproc*sim.getN_holes_collected() << " total holes collected out of " << transit_attempts_total << " total attempts.\n";
                }
                analysisfile << "Overall time-of-flight charge transport test results:\n";
                analysisfile << "Transit time is " << vector_avg(transit_times) << " ± " << vector_stdev(transit_times) << " s.\n";
                analysisfile << "Charge carrier mobility is " << vector_avg(mobilities) << " ± " << vector_stdev(mobilities) << " cm^2 V^-1 s^-1.\n";
            }
        }
        if(params_opv.Enable_dynamics_test){
            vector<double> times = sim.getDynamicsTransientTimes();
            vector<int> excitons_total = MPI_calculateVectorSum(sim.getDynamicsTransientExcitons());
            vector<int> electrons_total = MPI_calculateVectorSum(sim.getDynamicsTransientElectrons());
            vector<int> holes_total = MPI_calculateVectorSum(sim.getDynamicsTransientHoles());
            if(procid==0){
                ofstream transientfile;
                ss << "dynamics_average_transients.txt";
                transientfile.open(ss.str().c_str());
                ss.str("");
                transientfile << "Time (s),Exciton Density (cm^-3),Electron Density (cm^-3),Hole Density (cm^-3)" << endl;
                double volume_total = nproc*params_opv.Length*params_opv.Width*params_opv.Height*intpow(1e-7*params_opv.Unit_size,3);
                for(int i=0;i<(int)times.size();i++){
                    transientfile << times[i] << "," << excitons_total[i]/volume_total << "," << electrons_total[i]/volume_total << "," << holes_total[i]/volume_total << endl;
                }
                transientfile.close();
            }
        }
        if(params_opv.Enable_dynamics_test || params_opv.Enable_IQE_test){
            int excitons_created = sim.getN_excitons_created();
            int excitons_created_total;
            MPI_Reduce(&excitons_created,&excitons_created_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int excitons_created_donor = sim.getN_excitons_created((short)1);
            int excitons_created_donor_total;
            MPI_Reduce(&excitons_created_donor,&excitons_created_donor_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int excitons_created_acceptor = sim.getN_excitons_created((short)2);
            int excitons_created_acceptor_total;
            MPI_Reduce(&excitons_created_acceptor,&excitons_created_acceptor_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int excitons_dissociated = sim.getN_excitons_dissociated();
            int excitons_dissociated_total;
            MPI_Reduce(&excitons_dissociated,&excitons_dissociated_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int geminate_recombinations = sim.getN_geminate_recombinations();
            int geminate_recombinations_total;
            MPI_Reduce(&geminate_recombinations,&geminate_recombinations_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int bimolecular_recombinations = sim.getN_bimolecular_recombinations();
            int bimolecular_recombinations_total;
            MPI_Reduce(&bimolecular_recombinations,&bimolecular_recombinations_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int electrons_collected = sim.getN_electrons_collected();
            int electrons_collected_total;
            MPI_Reduce(&electrons_collected,&electrons_collected_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            int holes_collected = sim.getN_holes_collected();
            int holes_collected_total;
            MPI_Reduce(&holes_collected,&holes_collected_total,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            if(procid==0 && params_opv.Enable_dynamics_test){
                analysisfile << "Overall dynamics test results:\n";
            }
            if(procid==0 && params_opv.Enable_IQE_test){
                analysisfile << "Overall internal quantum efficiency test results:\n";
            }
            if(procid==0){
                analysisfile << excitons_created_total << " total excitons have been created.\n";
                analysisfile << excitons_created_donor_total << " excitons were created on donor sites.\n";
                analysisfile << excitons_created_acceptor_total << " excitons were created on acceptor sites.\n";
                analysisfile << 100*(double)excitons_dissociated_total/(double)excitons_created_total << "% of total excitons have dissociated.\n";
                analysisfile << 100*(double)geminate_recombinations_total/(double)excitons_dissociated_total << "% of total photogenerated charges were lost to geminate recombination.\n";
                analysisfile << 100*(double)bimolecular_recombinations_total/(double)excitons_dissociated_total << "% of total photogenerated charges were lost to bimolecular recombination.\n";
                analysisfile << 100*(double)(electrons_collected_total+holes_collected_total)/(2*(double)excitons_dissociated_total) << "% of total photogenerated charges were extracted.\n";
            }
            if(procid==0 && params_opv.Enable_IQE_test){
                analysisfile << "IQE = " << 100*(double)(electrons_collected_total+holes_collected_total)/(2*(double)excitons_created_total) << "%." << endl;
            }
        }
        if(procid==0){
            analysisfile.close();
        }
    }
    if(params_main.Enable_mpi){
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }
    return 0;
}

bool importParameters(ifstream* inputfile,Parameters_main& params_main,Parameters_OPV& params){
    string line;
    string var;
    size_t pos;
    vector<string> stringvars;
    bool error_status = false;
    while(inputfile->good()){
        getline(*inputfile,line);
        if((line.substr(0,2)).compare("--")!=0 && (line.substr(0,2)).compare("##")!=0){
            pos = line.find("/",0);
            var = line.substr(0,pos-1);
            stringvars.push_back(var);
        }
    }
    int i = 0;
    // General Parameters
    //enable_mpi
    params_main.Enable_mpi = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting mpi options" << endl;
        return false;
    }
    i++;
    //enable_logging
    params.Enable_logging = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting logging options" << endl;
        return false;
    }
    i++;
    //enable_periodic_x
    params.Enable_periodic_x = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting x-periodic boundary options" << endl;
        return false;
    }
    i++;
    //enable_periodic_y
    params.Enable_periodic_y = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting y-periodic boundary options" << endl;
        return false;
    }
    i++;
    //enable_periodic_z
    params.Enable_periodic_z = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting z-periodic boundary options" << endl;
        return false;
    }
    i++;
    params.Length = atoi(stringvars[i].c_str());
    i++;
    params.Width = atoi(stringvars[i].c_str());
    i++;
    params.Height = atoi(stringvars[i].c_str());
    i++;
    params.Unit_size = atof(stringvars[i].c_str());
    i++;
    params.Temperature = atoi(stringvars[i].c_str());
    i++;
    params.Recalc_cutoff = atoi(stringvars[i].c_str());
    i++;
    // Additional General Parameters
    params.Bias = atof(stringvars[i].c_str());
    i++;
    // Film Architecture Parameters
    params.Enable_neat = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling neat film architecture." << endl;
        return false;
    }
    i++;
    params.Enable_bilayer = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling bilayer film architecture." << endl;
        return false;
    }
    i++;
    params.Thickness_donor = atoi(stringvars[i].c_str());
    i++;
    params.Thickness_acceptor = atoi(stringvars[i].c_str());
    i++;
    params.Enable_random_blend = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling random blend film architecture." << endl;
        return false;
    }
    i++;
    params.Acceptor_conc = atof(stringvars[i].c_str());;
    i++;
    params_main.Enable_import_morphology_single = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling morphology import." << endl;
        return false;
    }
    i++;
    params_main.Morphology_filename = stringvars[i];
    i++;
    params_main.Enable_import_morphology_set = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling morphology set import." << endl;
        return false;
    }
    i++;
    params_main.Morphology_set_format = stringvars[i];
    i++;
    params_main.N_test_morphologies = atoi(stringvars[i].c_str());
    i++;
    params_main.N_morphology_set_size = atoi(stringvars[i].c_str());
    i++;
    // Test Parameters
    params.N_tests = atoi(stringvars[i].c_str());
    i++;
    params.Enable_exciton_diffusion_test = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling the exciton diffusion test." << endl;
        return false;
    }
    i++;
    params.Enable_ToF_test = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling the time-of-flight polaron transport test." << endl;
        return false;
    }
    i++;
    if(stringvars[i].compare("electron")==0){
        params.ToF_polaron_type = false;
    }
    else if(stringvars[i].compare("hole")==0){
        params.ToF_polaron_type = true;
    }
    else{
        cout << "Error setting polaron type for the time-of-flight test." << endl;
        return false;
    }
    i++;
    params.ToF_initial_polarons = atoi(stringvars[i].c_str());
    i++;
    params.ToF_transient_start = atof(stringvars[i].c_str());
    i++;
    params.ToF_transient_end = atof(stringvars[i].c_str());
    i++;
    params.ToF_pnts_per_decade = atoi(stringvars[i].c_str());
    i++;
    params.Enable_IQE_test = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling the internal quantum efficiency test." << endl;
        return false;
    }
    i++;
    params.IQE_time_cutoff = atof(stringvars[i].c_str());
    i++;
    params.Enable_dynamics_test = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error enabling the dynamics test." << endl;
        return false;
    }
    i++;
    params.Enable_dynamics_extraction = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting dynamics test extraction option." << endl;
        return false;
    }
    i++;
    params.Dynamics_initial_exciton_conc = atof(stringvars[i].c_str());
    i++;
    params.Dynamics_transient_start = atof(stringvars[i].c_str());
    i++;
    params.Dynamics_transient_end = atof(stringvars[i].c_str());
    i++;
    params.Dynamics_pnts_per_decade = atoi(stringvars[i].c_str());
    i++;
    // Exciton Parameters
    params.Exciton_generation_rate_donor = atof(stringvars[i].c_str());
    i++;
    params.Exciton_generation_rate_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Singlet_lifetime_donor = atof(stringvars[i].c_str());
    i++;
    params.Singlet_lifetime_acceptor = atof(stringvars[i].c_str());
    i++;
	params.Triplet_lifetime_donor = atof(stringvars[i].c_str());
	i++;
	params.Triplet_lifetime_acceptor = atof(stringvars[i].c_str());
	i++;
    params.R_singlet_hopping_donor = atof(stringvars[i].c_str());
    i++;
    params.R_singlet_hopping_acceptor = atof(stringvars[i].c_str());
    i++;
	params.R_triplet_hopping_donor = atof(stringvars[i].c_str());
	i++;
	params.R_triplet_hopping_acceptor = atof(stringvars[i].c_str());
	i++;
	params.Triplet_localization_donor = atof(stringvars[i].c_str());
	i++;
	params.Triplet_localization_acceptor = atof(stringvars[i].c_str());
	i++;
	params.Enable_FRET_triplet_annihilation = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting FRET triplet annihilation option." << endl;
		return false;
	}
	i++;
	params.R_exciton_exciton_annihilation_donor = atof(stringvars[i].c_str());
	i++;
	params.R_exciton_exciton_annihilation_acceptor = atof(stringvars[i].c_str());
	i++;
	params.R_exciton_polaron_annihilation_donor = atof(stringvars[i].c_str());
	i++;
	params.R_exciton_polaron_annihilation_acceptor = atof(stringvars[i].c_str());
	i++;
    params.FRET_cutoff = atoi(stringvars[i].c_str());
    i++;
    params.E_exciton_binding_donor = atof(stringvars[i].c_str());
    i++;
    params.E_exciton_binding_acceptor = atof(stringvars[i].c_str());
    i++;
    params.R_exciton_dissociation_donor = atof(stringvars[i].c_str());
    i++;
    params.R_exciton_dissociation_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Exciton_dissociation_cutoff = atoi(stringvars[i].c_str());
    i++;
	params.R_exciton_isc_donor = atof(stringvars[i].c_str());
	i++;
	params.R_exciton_isc_acceptor = atof(stringvars[i].c_str());
	i++;
	params.R_exciton_risc_donor = atof(stringvars[i].c_str());
	i++;
	params.R_exciton_risc_acceptor = atof(stringvars[i].c_str());
	i++;
	params.E_exciton_ST_donor = atof(stringvars[i].c_str());
	i++;
	params.E_exciton_ST_acceptor = atof(stringvars[i].c_str());
	i++;
    // Polaron Parameters
    params.Enable_phase_restriction = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting polaron phase restriction option." << endl;
        return false;
    }
    i++;
    params.R_polaron_hopping_donor = atof(stringvars[i].c_str());
    i++;
    params.R_polaron_hopping_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Polaron_localization_donor = atof(stringvars[i].c_str());
    i++;
    params.Polaron_localization_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Enable_miller_abrahams = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting Miller-Abrahams polaron hopping model options" << endl;
        return false;
    }
    i++;
    params.Enable_marcus = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting Marcus polaron hopping model options" << endl;
        return false;
    }
    i++;
    params.Reorganization_donor = atof(stringvars[i].c_str());
    i++;
    params.Reorganization_acceptor = atof(stringvars[i].c_str());
    i++;
    params.R_polaron_recombination = atof(stringvars[i].c_str());
    i++;
    params.Polaron_hopping_cutoff = atoi(stringvars[i].c_str());
    i++;
    params.Enable_gaussian_polaron_delocalization = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting Gaussian polaron delocalization option." << endl;
        return false;
    }
    i++;
    params.Polaron_delocalization_length = atof(stringvars[i].c_str());
    i++;
    // Lattice Parameters
    params.Homo_donor = atof(stringvars[i].c_str());
    i++;
    params.Lumo_donor = atof(stringvars[i].c_str());
    i++;
    params.Homo_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Lumo_acceptor = atof(stringvars[i].c_str());
    i++;
    //enable_gaussian_dos
    params.Enable_gaussian_dos = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting Gaussian DOS options" << endl;
        return false;
    }
    i++;
    params.Energy_stdev_donor = atof(stringvars[i].c_str());
    i++;
    params.Energy_stdev_acceptor = atof(stringvars[i].c_str());
    i++;
    //enable_exponential_dos
    params.Enable_exponential_dos = importBooleanParam(stringvars[i],error_status);
    if(error_status){
        cout << "Error setting Exponential DOS options" << endl;
        return false;
    }
    i++;
    params.Energy_urbach_donor = atof(stringvars[i].c_str());
    i++;
    params.Energy_urbach_acceptor = atof(stringvars[i].c_str());
    i++;
	//enable_correlated_disorder
	params.Enable_correlated_disorder = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting Correlated Disorder options" << endl;
		return false;
	}
	i++;
	params.Disorder_correlation_length = atof(stringvars[i].c_str());
	i++;
	//enable_gaussian_kernel
	params.Enable_gaussian_kernel = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting Correlated Disorder gaussian kernel options" << endl;
		return false;
	}
	i++;
	//enable_power_kernel
	params.Enable_power_kernel = importBooleanParam(stringvars[i], error_status);
	if (error_status) {
		cout << "Error setting Correlated Disorder gaussian kernel options" << endl;
		return false;
	}
	i++;
	params.Power_kernel_exponent = atoi(stringvars[i].c_str());
	i++;
    // Coulomb Calculation Parameters
    params.Dielectric_donor = atof(stringvars[i].c_str());
    i++;
    params.Dielectric_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Coulomb_cutoff = atoi(stringvars[i].c_str());
    i++;
    // Error checking
    if(params_main.Enable_import_morphology_set && !params_main.Enable_mpi){
        cout << "Error! MPI must be enabled in order to import a morphology set." << endl;
        return false;
    }
    return true;
}


