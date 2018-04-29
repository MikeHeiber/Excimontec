FLAGS = -Wall -Wextra -O3 -std=c++11
OBJS = main.o OSC_Sim.o Exciton.o Polaron.o Event.o Lattice.o Object.o Simulation.o Site.o Utils.o

Excimontec.exe : $(OBJS)
	ifeq ($(MPI),mpich2)
		mpicxx $(FLAGS) $(OBJS) -show -o Excimontec.exe
	endif
	ifeq ($(MPI),openmpi)
		mpicxx $(FLAGS) $(OBJS) -showme -o Excimontec.exe
	endif
	mpicxx $(FLAGS) $(OBJS) -o Excimontec.exe

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
	
clean:
	\rm *.o *~ Excimontec.exe