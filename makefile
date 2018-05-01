# Copyright (c) 2018 Michael C. Heiber
# This source file is part of the Excimontec project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

ifeq ($(lastword $(subst /, ,$(CXX))),g++)
	FLAGS = -Wall -Wextra -O3 -std=c++11 -I. -Isrc
endif
ifeq ($(lastword $(subst /, ,$(CXX))),pgc++)
	FLAGS = -O2 -fastsse -Mvect -std=c++11 -Mdalign -Munroll -Mipa=fast -Kieee -m64 -I. -Isrc
endif

OBJS = src/OSC_Sim.o src/Exciton.o src/Polaron.o src/Event.o src/Lattice.o src/Object.o src/Simulation.o src/Site.o src/Utils.o

all : Excimontec.exe
ifndef FLAGS
	$(error Valid compiler not detected.)
endif

Excimontec.exe : src/main.o $(OBJS)
	mpicxx -v $(FLAGS) src/main.o $(OBJS) -o Excimontec.exe

src/main.o : src/main.cpp src/OSC_Sim.h src/Exciton.h src/Polaron.h KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c src/main.cpp -o $@
	
src/OSC_Sim.o : src/OSC_Sim.h src/OSC_Sim.cpp src/Exciton.h src/Polaron.h KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c src/OSC_Sim.cpp -o $@

src/Exciton.o : src/Exciton.h src/Exciton.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c src/Exciton.cpp -o $@

src/Polaron.o : src/Polaron.h src/Polaron.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c src/Polaron.cpp -o $@

src/Event.o : KMC_Lattice/Event.h KMC_Lattice/Event.cpp KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Simulation.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Event.cpp -o $@

src/Lattice.o : KMC_Lattice/Lattice.h KMC_Lattice/Lattice.cpp KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Lattice.cpp -o $@

src/Object.o : KMC_Lattice/Object.h KMC_Lattice/Object.cpp KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Object.cpp -o $@

src/Simulation.o : KMC_Lattice/Simulation.h KMC_Lattice/Simulation.cpp KMC_Lattice/Event.h KMC_Lattice/Lattice.h KMC_Lattice/Object.h KMC_Lattice/Site.h KMC_Lattice/Utils.h
	mpicxx $(FLAGS) -c KMC_Lattice/Simulation.cpp -o $@
	
src/Site.o : KMC_Lattice/Site.h KMC_Lattice/Site.cpp
	mpicxx $(FLAGS) -c KMC_Lattice/Site.cpp -o $@
	
src/Utils.o : KMC_Lattice/Utils.h KMC_Lattice/Utils.cpp
	mpicxx $(FLAGS) -c KMC_Lattice/Utils.cpp -o $@

#
# Testing Section using googletest
#

GTEST_DIR = googletest/googletest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

ifeq ($(lastword $(subst /, ,$(CXX))),g++)
	GTEST_FLAGS = -isystem $(GTEST_DIR)/include -pthread
endif
ifeq ($(lastword $(subst /, ,$(CXX))),pgc++)
	GTEST_FLAGS = -I$(GTEST_DIR)/include
endif
ifndef FLAGS
	$(error Valid compiler not detected.)
endif

testing/gtest-all.o : $(GTEST_SRCS_)
	mpicxx $(GTEST_FLAGS) -I$(GTEST_DIR) $(FLAGS) -c $(GTEST_DIR)/src/gtest-all.cc -o $@
			
testing/test.o : testing/test.cpp $(GTEST_HEADERS) src/OSC_Sim.h src/Exciton.h KMC_Lattice/Utils.h
	mpicxx $(GTEST_FLAGS) $(FLAGS) -c testing/test.cpp -o $@

testing/Excimontec_tests.exe : testing/test.o testing/gtest-all.o $(OBJS)
	mpicxx $(GTEST_FLAGS) $(FLAGS) -lpthread $^ -o $@
			
test : testing/Excimontec_tests.exe
ifndef FLAGS
	$(error Valid compiler not detected.)
endif
	
clean:
	\rm src/*.o *~ Excimontec.exe testing/*.o testing/Excimontec_tests.exe
