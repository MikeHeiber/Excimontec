#!/bin/bash
#SBATCH -J Excimontec # Job name
#SBATCH -p partition_name
#SBATCH -n 48 # Number of tasks
#SBATCH -t 01:00:00 # Maximum walltime
#SBATCH --cpus-per-task=1

version_num=v1.0-beta.3
ParameterNum=default

# Setup job directory
mkdir $SLURM_JOB_ID
cd $SLURM_JOB_ID
cp ../Excimontec.exe ./Excimontec.exe
cp ../parameters_$ParameterNum.txt ./parameters_$ParameterNum.txt

# Execute Excimontec Simulation
mpiexec -n 48 Excimontec.exe parameters_$ParameterNum.txt > output.txt

# Cleanup
rm -f Excimontec.exe
tar -zcf $SLURM_JOB_ID.tar.gz $SLURM_JOB_ID
cp $SLURM_JOB_ID.tar.gz ../$SLURM_JOB_ID.tar.gz
rm -f $SLURM_JOB_ID.tar.gz