# Copyright (c) 2018 Michael C. Heiber
# This source file is part of the Excimontec project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec
# makefile for GNU make

FLAGS = -Wall -Wextra -O3 -std=c++11
OBJS = main.o OSC_Sim.o Exciton.o Polaron.o Event.o Lattice.o Object.o Simulation.o Site.o Utils.o

Excimontec.exe : $(OBJS)
	mpicxx -v $(FLAGS) $(OBJS) -o Excimontec.exe

main.o : main.cpp OSC_Sim.h Exciton.h Polaron.h KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c main.cpp
	
OSC_Sim.o : OSC_Sim.h OSC_Sim.cpp Exciton.h Polaron.h KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c OSC_Sim.cpp

Exciton.o : Exciton.h Exciton.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c Exciton.cpp

Polaron.o : Polaron.h Polaron.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c Polaron.cpp

Event.o : KMC_Lattice/Event.h KMC_Lattice/Event.cpp KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Event.cpp

Lattice.o : KMC_Lattice/Lattice.h KMC_Lattice/Lattice.cpp KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Lattice.cpp

Object.o : KMC_Lattice/Object.h KMC_Lattice/Object.cpp KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Object.cpp

Simulation.o : KMC_Lattice/Simulation.h KMC_Lattice/Simulation.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Simulation.cpp
	
Site.o : KMC_Lattice/Site.h KMC_Lattice/Site.cpp
	mpicxx $(FLAGS) -c KMC_Lattice/Site.cpp
	
Utils.o : KMC_Lattice/Utils.h KMC_Lattice/Utils.cpp
	mpicxx $(FLAGS) -c KMC_Lattice/Utils.cpp

# Testing Section using googletest
GTEST_DIR = googletest/googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
GTEST_FLAGS = -isystem $(GTEST_DIR)/include -pthread

gtest-all.o : $(GTEST_SRCS_)
	mpicxx $(GTEST_FLAGS) -I$(GTEST_DIR) $(FLAGS) -c $(GTEST_DIR)/src/gtest-all.cc
			
test.o : testing/test.cpp $(GTEST_HEADERS) OSC_Sim.h Exciton.h KMC_Lattice/Utils.h
	mpicxx $(GTEST_FLAGS) $(FLAGS) -c testing/test.cpp

Excimontec_tests.exe : test.o gtest-all.o
	mpicxx $(GTEST_FLAGS) $(FLAGS) -lpthread $^ -o $@
			
# Run the tests
test : Excimontec_tests.exe

	
clean:
	\rm *.o *~ Excimontec.exe Excimontec_tests.exe