#!/bin/bash
#SBATCH -J Excimontec # Job name
#SBATCH -n 16 # Number of tasks
#SBATCH -t 00:30:00 # Maximum walltime
#SBATCH --cpus-per-task=1

# Setup job directory
mkdir $SLURM_JOB_ID
cd $SLURM_JOB_ID
cp ../../Excimontec.exe ./Excimontec.exe
cp ../parameters_test1.txt ./parameters_test1.txt

# Execute Excimontec Simulation
mpiexec -n 16 Excimontec.exe parameters_test1.txt > test1_output.txt

# Cleanup
cp test1_output.txt ../test1_output.txt
rm -f Excimontec.exe
