# Excimontec
The goal of this project is to develop an open-source lattice KMC simulation software package for modeling organic semiconductor materials and devices, such as OPVs, OLEDs, and more.  The software is being developed in C++11 and optimized for efficient execution on high performance computing clusters using MPI.  This software package uses object-oriented design and extends the [KMC_Lattice](https://github.com/MikeHeiber/KMC_Lattice) framework.

## Current Status:
The current version (Excimontec v0.1-alpha) is built with KMC_Lattice v1.1-alpha.1 and allows the user to perform several simulation tests relevant for OPV and OLED devices.  NOTE: This software is still in the alpha phase of development, and as such, there may still be bugs that need to be squashed.  Please report any bugs or submit feature requests in the Issues section.

#### Implemented Major Features:
- Choose between several film architectures, including a neat film, bilayer film, or random blend film.
- Donor and acceptor materials can take on an uncorrelated Gaussian DOS or an exponential DOS model.
- Time-of-flight charge transport simulations of electrons or holes can be performed on neat or random blend films.
- Exciton diffusion simulations can be performed on any film architecture.
- Internal quantum efficiency simulations can be performed on bilayer or random blend films.
- Adjustable periodic boundary conditions in all three directions allow users to perform 1D, 2D, or 3D simulations.
- Choose between Miller-Abrahams or Marcus models for polaron hopping.

#### Features coming soon:
- Import model bulk heterojunction morpholgies generated with the [Ising_OPV](https://github.com/MikeHeiber/Ising_OPV) software tool.
- Bimolecular charge recombination test

## How to try Excimontec?
#### Compiling
Compiling requires an MPI library for parallel processing.

More information about these packages can be found here:
- http://www.mpich.org/, http://www.open-mpi.org/, http://mvapich.cse.ohio-state.edu/

#### Usage
Excimontec.exe takes one required input argument, which is the filename of the input parameter file.

An example parameter file is provided with parameters_default.txt

As an example, to create a single simulation instance on a single processor, the command is:
>   Excimontec.exe parameters_default.txt

To run in a parallel processing environment and create 10 simulations on 10 processors to gather more statistics, an example run command is:
>    mpiexec -n 10 Excimontec.exe parameters_default.txt

MPI execution commands can be implemented into batch scripts for running Excimontec in a supercomputing environment.

#### Output
Excimontec will create several output files depending which test is chosen:
- results#.txt -- This text file will contain the results for each processor where the # will be replaced by the processor ID.
- analysis_summary.txt -- When MPI is enabled, this text file will contain average final results from all of the processors.
- ToF_average_transients.txt -- When performing a time-of-flight charge transport test, calculated current transients, mobility relaxation transients, and energy relaxation transients will be output to this file.
- ToF_transit_time_dist.txt -- When performing a time-of-flight charge transport test, the resulting polaron transit time probability distribution will be output to this file.
