#!/bin/bash

# Unit tests
./Excimontec_tests.exe
./Version_tests.exe

# System tests
sbatch system_test1.sh
