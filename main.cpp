// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#include "KMC_Lattice/Utils.h"
#include "OSC_Sim.h"
#include "mpi.h"
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

using namespace std;

struct Parameters_main{
    bool Enable_mpi;
};

//Declare Functions
bool importParameters(ifstream * inputfile,Parameters_main& params_main,Parameters_OPV& params);

int main(int argc,char *argv[]){
    string version = "v1.0";
    // Parameters
    bool End_sim = false;
    // File declaration
    ifstream parameterfile;
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
    if(!parameterfile){
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
    if(params_main.Enable_mpi){
        cout << "Initializing MPI options... ";
        MPI::Init(argc,argv);
        nproc = MPI::COMM_WORLD.Get_size();
        procid = MPI::COMM_WORLD.Get_rank();
        cout << "MPI initialization complete!" << endl;
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
    OSC_Sim sim(params_opv,procid);
    cout << procid << ": Simulation initialization complete" << endl;
    // Begin Simulation loop
    while(!End_sim){
        success = sim.executeNextEvent();
        if(!success){
            cout << procid << ": Event execution failed, simulation will now terminate." << endl;
            break;
        }
        // Check if simulation has finished
        End_sim = sim.checkFinished();
        // Output status
        if(sim.getN_events_executed()%10000==0){
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
    elapsedtime = difftime(time_end,time_start);
    // Output simulation results for each processor
    ss << "results" << procid << ".txt";
    resultsfile.open(ss.str().c_str());
    ss.str("");
    resultsfile << "KMC_Lattice_example " << version << " Results:\n";
    resultsfile << "Calculation time elapsed is " << (double)elapsedtime/60 << " minutes.\n";
    resultsfile << sim.getTime() << " seconds have been simulated.\n";
    resultsfile << sim.getN_events_executed() << " events have been executed.\n";
    resultsfile << sim.getN_excitons_created() << " excitons have been created.\n";
    if(params_opv.Enable_exciton_diffusion_test){
        resultsfile << "Exciton diffusion test results:\n";
        resultsfile << "Exciton Diffusion Length is " << sim.calculateDiffusionLength_avg() << " ± " << sim.calculateDiffusionLength_stdev() << " nm\n";
    }
    resultsfile << endl;
    resultsfile.close();
    // Output overall analysis results from all processors
    if(params_main.Enable_mpi){
        vector<double> diffusion_data;
        if(params_opv.Enable_exciton_diffusion_test){
            diffusion_data = calculateAverageVector(sim.getDiffusionData(),procid,nproc);
        }
        if(procid==0){
            ss << "analysis_summary.txt";
            analysisfile.open(ss.str().c_str());
            ss.str("");
            analysisfile << "KMC_Lattice_example " << version << " Results Summary:" << endl;
            analysisfile << nproc*sim.getN_excitons_recombined() << " total excitons tested." << endl;
            if(params_opv.Enable_exciton_diffusion_test){
                analysisfile << "Overall exciton diffusion test results:\n";
                analysisfile << "Exciton diffusion length is " << vector_avg(diffusion_data) << " ± " << vector_stdev(diffusion_data) << " nm\n";
            }
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
    if(stringvars[i].compare("true")==0){
        params_main.Enable_mpi = true;
    }
    else if(stringvars[i].compare("false")==0){
        params_main.Enable_mpi = false;
    }
    else{
        cout << "Error setting mpi options" << endl;
        return false;
    }
    i++;
    //enable_logging
    if(stringvars[i].compare("true")==0){
        params.Enable_logging = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_logging = false;
    }
    else{
        cout << "Error setting logging options" << endl;
        return false;
    }
    i++;
    //enable_periodic_x
    if(stringvars[i].compare("true")==0){
        params.Enable_periodic_x = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_periodic_x = false;
    }
    else{
        cout << "Error setting x-periodic boundary options" << endl;
        return false;
    }
    i++;
    //enable_periodic_y
    if(stringvars[i].compare("true")==0){
        params.Enable_periodic_y = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_periodic_y = false;
    }
    else{
        cout << "Error setting y-periodic boundary options" << endl;
        return false;
    }
    i++;
    //enable_periodic_z
    if(stringvars[i].compare("true")==0){
        params.Enable_periodic_z = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_periodic_z = false;
    }
    else{
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
    if(stringvars[i].compare("true")==0){
        params.Enable_neat = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_neat = false;
    }
    else{
        cout << "Error enabling neat film architecture." << endl;
        return false;
    }
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_bilayer = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_bilayer = false;
    }
    else{
        cout << "Error enabling bilayer film architecture." << endl;
        return false;
    }
    i++;
    params.Thickness_donor = atoi(stringvars[i].c_str());
    i++;
    params.Thickness_acceptor = atoi(stringvars[i].c_str());
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_random_blend = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_random_blend = false;
    }
    else{
        cout << "Error enabling random blend film architecture." << endl;
        return false;
    }
    i++;
    params.Acceptor_conc = atof(stringvars[i].c_str());;
    i++;
    // Test Parameters
    params.N_tests = atoi(stringvars[i].c_str());
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_exciton_diffusion_test = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_exciton_diffusion_test = false;
    }
    else{
        cout << "Error enabling the exciton diffusion test." << endl;
        return false;
    }
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_ToF_test = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_ToF_test = false;
    }
    else{
        cout << "Error enabling the time-of-flight polaron transport test." << endl;
        return false;
    }
    i++;
    params.ToF_initial_polarons = atoi(stringvars[i].c_str());
    i++;
    params.ToF_start_time = atof(stringvars[i].c_str());
    i++;
    params.ToF_end_time = atof(stringvars[i].c_str());
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_IQE_test = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_IQE_test = false;
    }
    else{
        cout << "Error enabling the internal quantum efficiency test." << endl;
        return false;
    }
    i++;
    // Exciton Parameters
    params.Exciton_generation_rate_donor = atof(stringvars[i].c_str());
    i++;
    params.Exciton_generation_rate_acceptor = atof(stringvars[i].c_str());
    i++;
    params.Exciton_lifetime_donor = atof(stringvars[i].c_str());
    i++;
    params.Exciton_lifetime_acceptor = atof(stringvars[i].c_str());
    i++;
    params.R_exciton_hopping_donor = atof(stringvars[i].c_str());
    i++;
    params.R_exciton_hopping_acceptor = atof(stringvars[i].c_str());
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
    // Polaron Parameters
    if(stringvars[i].compare("true")==0){
        params.Enable_phase_restriction = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_phase_restriction = false;
    }
    else{
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
    params.Polaron_hopping_cutoff = atoi(stringvars[i].c_str());
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_miller_abrahams = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_miller_abrahams = false;
    }
    else{
        cout << "Error setting Miller-Abrahams polaron hopping model options" << endl;
        return false;
    }
    i++;
    if(stringvars[i].compare("true")==0){
        params.Enable_marcus = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_marcus = false;
    }
    else{
        cout << "Error setting Marcus polaron hopping model options" << endl;
        return false;
    }
    i++;
    params.Reorganization_energy_donor = atof(stringvars[i].c_str());
    i++;
    params.Reorganization_energy_acceptor = atof(stringvars[i].c_str());
    i++;
    // Energetic Disorder Parameters
    //enable_gaussian_dos
    if(stringvars[i].compare("true")==0){
        params.Enable_gaussian_dos = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_gaussian_dos = false;
    }
    else{
        cout << "Error setting Gaussian DOS options" << endl;
        return false;
    }
    i++;
    params.Energy_stdev_donor = atof(stringvars[i].c_str());
    i++;
    params.Energy_stdev_acceptor = atof(stringvars[i].c_str());
    i++;
    //enable_exponential_dos
    if(stringvars[i].compare("true")==0){
        params.Enable_exponential_dos = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_exponential_dos = false;
    }
    else{
        cout << "Error setting Exponential DOS options" << endl;
        return false;
    }
    i++;
    params.Energy_urbach_donor = atof(stringvars[i].c_str());
    i++;
    params.Energy_urbach_acceptor = atof(stringvars[i].c_str());
    i++;
    // Coulomb Calculation Parameters
    params.Coulomb_cutoff = atoi(stringvars[i].c_str());
    i++;
    params.Dielectric_donor = atof(stringvars[i].c_str());
    i++;
    params.Dielectric_acceptor = atof(stringvars[i].c_str());
    i++;
    // Error checking
    if(!params.Length>0 || !params.Width>0 || !params.Height>0){
        cout << "Error! All lattice dimensions must be greater than zero." << endl;
        return false;
    }
    if(!params.Unit_size>0){
        cout << "Error! The lattice unit size must be greater than zero." << endl;
        return false;
    }
    if(!params.Temperature>0){
        cout << "Error! The temperature must be greater than zero." << endl;
        return false;
    }
    if(params.Recalc_cutoff<params.FRET_cutoff){
        cout << "Error! The event recalculation cutoff radius must not be less than the FRET cutoff radius." << endl;
        return false;
    }
    if(params.Recalc_cutoff<params.Polaron_hopping_cutoff){
        cout << "Error! The event recalculation cutoff radius must not be less than the polaron hopping cutoff radius." << endl;
        return false;
    }
    if(!params.N_tests>0){
        cout << "Error! The number of tests must be greater than zero." << endl;
        return false;
    }
    if(params.Enable_gaussian_dos && params.Enable_exponential_dos){
        cout << "Error! The Gaussian and exponential disorder models cannot both be enabled." << endl;
        return false;
    }
    if(params.Enable_gaussian_dos && (params.Energy_stdev_donor<0 || params.Energy_stdev_acceptor<0)){
        cout << "Error! When using the Gaussian disorder model, the standard deviation cannot be negative." << endl;
        return false;
    }
    if(params.Enable_exponential_dos && (params.Energy_urbach_donor<0 || params.Energy_urbach_acceptor<0)){
        cout << "Error! When using the exponential disorder model, the Urbach energy cannot be negative." << endl;
        return false;
    }
    return true;
}


