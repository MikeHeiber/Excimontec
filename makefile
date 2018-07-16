# Copyright (c) 2018 Michael C. Heiber
# This source file is part of the Excimontec project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

ifeq ($(lastword $(subst /, ,$(CXX))),g++)
	FLAGS += -Wall -Wextra -O3 -std=c++11 
	ifdef ENABLE_OMP
		FLAGS += -openmp
	endif
	FLAGS += -I. -Isrc -IKMC_Lattice/src
endif
ifeq ($(lastword $(subst /, ,$(CXX))),pgc++)
	FLAGS += -O2 -fastsse -Mvect -std=c++11 -Mdalign -Munroll -Mipa=fast -Kieee -m64 
	ifdef ENABLE_OMP
		FLAGS += -mp
	endif
	FLAGS += -I. -Isrc -IKMC_Lattice/src
endif


OBJS = src/OSC_Sim.o src/Exciton.o src/Polaron.o

all : Excimontec.exe
ifndef FLAGS
	$(error Valid compiler not detected.)
endif

Excimontec.exe : src/main.o $(OBJS) KMC_Lattice/libKMC.a
	mpicxx $(FLAGS) $^ -o $@
	
KMC_Lattice/libKMC.a : KMC_Lattice/src/*.h
	$(MAKE) -C KMC_Lattice 

src/main.o : src/main.cpp src/OSC_Sim.h src/Exciton.h src/Polaron.h
	mpicxx $(FLAGS) -c $< -o $@
	
src/OSC_Sim.o : src/OSC_Sim.cpp src/OSC_Sim.h src/Exciton.h src/Polaron.h
	mpicxx $(FLAGS) -c $< -o $@

src/Exciton.o : src/Exciton.cpp src/Exciton.h
	mpicxx $(FLAGS) -c $< -o $@

src/Polaron.o : src/Polaron.cpp src/Polaron.h
	mpicxx $(FLAGS) -c $< -o $@

#
# Testing Section using googletest
#

ifndef FLAGS
	$(error Valid compiler not detected.)
endif
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

test_coverage : FLAGS = -fprofile-arcs -ftest-coverage -std=c++11 -Wall -Wextra -I. -Isrc -IKMC_Lattice/src
test_coverage : test/Excimontec_tests.exe

test : test/Excimontec_tests.exe	
	
test/Excimontec_tests.exe : test/test.o test/gtest-all.o $(OBJS) KMC_Lattice/libKMC.a
	mpicxx $(GTEST_FLAGS) $(FLAGS) $^ KMC_Lattice/libKMC.a -lpthread -o $@

test/gtest-all.o : $(GTEST_SRCS_)
	mpicxx $(GTEST_FLAGS) -I$(GTEST_DIR) $(FLAGS) -c $(GTEST_DIR)/src/gtest-all.cc -o $@
			
test/test.o : test/test.cpp $(GTEST_HEADERS) $(OBJS)
	mpicxx $(GTEST_FLAGS) $(FLAGS) -c $< -o $@
	
clean:
	-rm src/*.o src/*.gcno* src/*.gcda test/*.o test/*.gcno* test/*.gcda *~ Excimontec.exe test/Excimontec_tests.exe
