// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef POLARON_H
#define POLARON_H

#include "KMC_Lattice/Utils.h"
#include "KMC_Lattice/Object.h"
#include "KMC_Lattice/Event.h"
#include <string>

class Polaron : public Object{
    public:
        static const std::string name;
        Polaron(const double time,const int tag_num,const Coords& start_coords,const bool polaron_charge) : Object(time,tag_num,start_coords){charge = polaron_charge;}
        bool getCharge() const{return charge;}
		std::string getName() const{return name;}
    private:
        bool charge; // false represents negative charge, true represents positive charge
};

class Polaron_Hop : public Event{
    public:
        static const std::string name;
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,Simulation* sim_ptr){
            // Calculates hopping using the Miller-Abrahams model
            double rate = prefactor*exp(-2.0*localization*distance);
            if(E_delta>0){
                rate *= exp(-E_delta/(Utils::K_b*sim_ptr->getTemp()));
            }
            Event::calculateExecutionTime(rate,sim_ptr);
        }
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,const double reorganization, Simulation* sim_ptr){
            // Calculates hopping using the Marcus model
            double rate = (prefactor/sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization+E_delta,2)/(4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
            Event::calculateExecutionTime(rate,sim_ptr);
        }
		std::string getName() const{return name;}
    private:

};

class Polaron_Recombination : public Event{
    public:
        static const std::string name;
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta, Simulation* sim_ptr){
            // Calculates recombination using the Miller-Abrahams model
            double rate = prefactor*exp(-2.0*localization*distance);
            if(E_delta>0){
                rate *= exp(-E_delta/ (Utils::K_b*sim_ptr->getTemp()));
            }
            Event::calculateExecutionTime(rate,sim_ptr);
        }
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,const double reorganization, Simulation* sim_ptr){
            // Calculates recombination using the Marcus model
            double rate = (prefactor/sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization+E_delta,2)/(4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
            Event::calculateExecutionTime(rate,sim_ptr);
        }
		std::string getName() const{return name;}
    private:

};

class Polaron_Extraction : public Event{
    public:
        static const std::string name;
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta, Simulation* sim_ptr){
            // Calculates extraction using the Miller-Abrahams model
            double rate = prefactor*exp(-2.0*localization*distance);
            if(E_delta>0){
                rate *= exp(-E_delta/(Utils::K_b*sim_ptr->getTemp()));
            }
            Event::calculateExecutionTime(rate,sim_ptr);
        }
		std::string getName() const{return name;}
    private:
};

#endif // POLARON_H
