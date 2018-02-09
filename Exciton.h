// Copyright (c) 2018 Michael C. Heiber
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
        static const std::string object_type;
        Exciton(const double time,const int tag_num,const Coords& start_coords) : Object(time,tag_num,start_coords){}
        void flipSpin(){spin_state = !spin_state;}
		std::string getObjectType() const{return object_type;}
        bool getSpin() const{return spin_state;}
        void setSpin(bool spin_state_new){spin_state = spin_state_new;}
    private:
        bool spin_state; // false represents triplet state, true represents singlet state
};

class Exciton_Creation : public Event{
    public:
        static const std::string event_type;
		Exciton_Creation() : Event() {}
		Exciton_Creation(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		std::string getEventType() const{return event_type;}
    private:
};

class Exciton_Hop : public Event{
    public:
        static const std::string event_type;
		Exciton_Hop() : Event() {}
		Exciton_Hop(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		// Singlet FRET hop
		void calculateExecutionTime(const double prefactor, const double distance, const double E_delta) {
			double rate = prefactor*Utils::intpow(1.0 / distance, 6);
			if (E_delta > 0) {
				rate *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
			Event::calculateExecutionTime(rate);
		}
		// Triplet Dexter hop
		void calculateExecutionTime(const double prefactor, const double localization, const double distance, const double E_delta) {
			double rate = prefactor*exp(-2.0*localization*distance);
			if (E_delta>0) {
				rate *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
			Event::calculateExecutionTime(rate);
		}
		std::string getEventType() const{return event_type;}
    private:
};

class Exciton_Recombination : public Event{
    public:
        static const std::string event_type;
		Exciton_Recombination() : Event() {}
		Exciton_Recombination(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		std::string getEventType() const{return event_type;}
    private:
};

class Exciton_Dissociation : public Event{
    public:
        static const std::string event_type;
		Exciton_Dissociation() : Event() {}
		Exciton_Dissociation(Simulation* simulation_ptr) : Event(simulation_ptr) {}
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta){
            // Calculates dissociation using the Miller-Abrahams model
            double rate = prefactor*exp(-2.0*localization*distance);
            if(E_delta>0){
                rate *= exp(-E_delta/(Utils::K_b*sim_ptr->getTemp()));
            }
            Event::calculateExecutionTime(rate);
        }
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,const double reorganization){
            // Calculates dissociation using the Marcus model
            double rate = (prefactor/sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization+E_delta,2)/(4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
            Event::calculateExecutionTime(rate);
        }
		std::string getEventType() const{return event_type;}
    private:
};

class Exciton_Intersystem_Crossing : public Event {
    public:
        static const std::string event_type;
		Exciton_Intersystem_Crossing() : Event() {}
		Exciton_Intersystem_Crossing(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		void calculateExecutionTime(const double prefactor, const double E_delta) {
			double rate = prefactor;
			if (E_delta>0) {
				rate *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
			Event::calculateExecutionTime(rate);
		}
		std::string getEventType() const{return event_type;}
    private:
};

class Exciton_Exciton_Annihilation : public Event{
    public:
		static const std::string event_type;
		Exciton_Exciton_Annihilation() : Event() {}
		Exciton_Exciton_Annihilation(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		void calculateExecutionTime(const double prefactor, const double distance) {
			Event::calculateExecutionTime(prefactor*Utils::intpow(1.0 / distance, 6));
		}
		void calculateExecutionTime(const double prefactor, const double localization, const double distance) {
			Event::calculateExecutionTime(prefactor*exp(-2.0*localization*distance));
		}
		std::string getEventType() const { return event_type; }
    private:
};

class Exciton_Polaron_Annihilation : public Event{
    public:
		static const std::string event_type;
		Exciton_Polaron_Annihilation() : Event() {}
		Exciton_Polaron_Annihilation(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		void calculateExecutionTime(const double prefactor, const double distance) {
			Event::calculateExecutionTime(prefactor*Utils::intpow(1.0 / distance, 6));
		}
		void calculateExecutionTime(const double prefactor, const double localization, const double distance) {
			Event::calculateExecutionTime(prefactor*exp(-2.0*localization*distance));
		}
		std::string getEventType() const { return event_type; }
    private:
};

#endif // EXCITON_H
