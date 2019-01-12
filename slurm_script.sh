#!/bin/bash
# Copyright (c) 2017-2019 Michael C. Heiber
# This source file is part of the Excimontec project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#SBATCH -J Excimontec # Job name
#SBATCH -p partition_name
#SBATCH -n 48 # Number of tasks
#SBATCH -t 01:00:00 # Maximum walltime
#SBATCH --cpus-per-task=1

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