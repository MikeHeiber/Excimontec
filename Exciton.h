// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCITON_H
#define EXCITON_H

#include "KMC_Lattice/Utils.h"
#include "KMC_Lattice/Object.h"
#include "KMC_Lattice/Event.h"
#include "KMC_Lattice/Simulation.h"
#include <string>

class Exciton : public Object{
    public:
        static const std::string name;
        Exciton(const double time,const int tag_num,const Coords& start_coords) : Object(time,tag_num,start_coords){}
        void flipSpin(){spin = !spin;}
		std::string getName() const{return name;}
        bool getSpin() const{return spin;}
        void setSpin(bool spin_state){spin = spin_state;}
    private:
        bool spin; // false represents triplet state, true represents singlet state
};

class Exciton_Creation : public Event{
    public:
        static const std::string name;
		std::string getName() const{return name;}
    private:
};

class Exciton_Hop : public Event{
    public:
        static const std::string name;
		// Singlet FRET hop
		void calculateExecutionTime(const double prefactor, const double distance, const double E_delta, Simulation* sim_ptr) {
			double rate = prefactor*Utils::intpow(1.0 / distance, 6);
			if (E_delta > 0) {
				rate *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
			Event::calculateExecutionTime(rate, sim_ptr);
		}
		// Triplet Dexter hop
		void calculateExecutionTime(const double prefactor, const double localization, const double distance, const double E_delta, Simulation* sim_ptr) {
			double rate = prefactor*exp(-2.0*localization*distance);
			if (E_delta>0) {
				rate *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
			Event::calculateExecutionTime(rate, sim_ptr);
		}
		std::string getName() const{return name;}
    private:
};

class Exciton_Recombination : public Event{
    public:
        static const std::string name;
		std::string getName() const{return name;}
    private:
};

class Exciton_Dissociation : public Event{
    public:
        static const std::string name;
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,Simulation* sim_ptr){
            // Calculates dissociation using the Miller-Abrahams model
            double rate = prefactor*exp(-2.0*localization*distance);
            if(E_delta>0){
                rate *= exp(-E_delta/(Utils::K_b*sim_ptr->getTemp()));
            }
            Event::calculateExecutionTime(rate,sim_ptr);
        }
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,const double reorganization,Simulation* sim_ptr){
            // Calculates dissociation using the Marcus model
            double rate = (prefactor/sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization+E_delta,2)/(4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
            Event::calculateExecutionTime(rate,sim_ptr);
        }
		std::string getName() const{return name;}
    private:
};

class Exciton_Intersystem_Crossing : public Event {
    public:
        static const std::string name;
		void calculateExecutionTime(const double prefactor, const double E_delta, Simulation* sim_ptr) {
			double rate = prefactor;
			if (E_delta>0) {
				rate *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
			Event::calculateExecutionTime(rate, sim_ptr);
		}
		std::string getName() const{return name;}
    private:
};

class Exciton_Exciton_Annihilation : public Event{
    public:
		static const std::string name;
		void calculateExecutionTime(const double prefactor, const double distance, Simulation* sim_ptr) {
			Event::calculateExecutionTime(prefactor*Utils::intpow(1.0 / distance, 6), sim_ptr);
		}
		void calculateExecutionTime(const double prefactor, const double localization, const double distance, Simulation* sim_ptr) {
			Event::calculateExecutionTime(prefactor*exp(-2.0*localization*distance), sim_ptr);
		}
		std::string getName() const { return name; }
    private:
};

class Exciton_Polaron_Annihilation : public Event{
    public:
		static const std::string name;
		void calculateExecutionTime(const double prefactor, const double distance, Simulation* sim_ptr) {
			Event::calculateExecutionTime(prefactor*Utils::intpow(1.0 / distance, 6), sim_ptr);
		}
		void calculateExecutionTime(const double prefactor, const double localization, const double distance, Simulation* sim_ptr) {
			Event::calculateExecutionTime(prefactor*exp(-2.0*localization*distance), sim_ptr);
		}
		std::string getName() const { return name; }
    private:
};

#endif // EXCITON_H
