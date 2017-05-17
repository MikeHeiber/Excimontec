// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCITON_H
#define EXCITON_H

#include "KMC_Lattice/Utils.h"
#include "KMC_Lattice/Object.h"
#include "KMC_Lattice/Event.h"
#include <string>

using namespace std;

class Exciton : public Object{
    public:
        static const string name;
        Exciton(const double time,const int tag_num,const Coords& start_coords) : Object(time,tag_num,start_coords){}
        void flipSpin(){spin = !spin;}
        string getName() const{return name;}
        bool getSpin() const{return spin;}
        void setSpin(bool spin_state){spin = spin_state;}
    private:
        bool spin; // false represents triplet state, true represents singlet state
};

class Exciton_Creation : public Event{
    public:
        static const string name;
        string getName() const{return name;}
    private:


};

class Exciton_Hop : public Event{
    using Event::calculateExecutionTime;
    public:
        static const string name;
        void calculateExecutionTime(const double prefactor,const double distance,const double E_delta,const double temp,const double current_time){
            double rate = prefactor*intpow(1/distance,6);
            if(E_delta>0){
                rate *= exp(-E_delta/(K_b*temp));
            }
            calculateExecutionTime(rate,current_time);
        }
        string getName() const{return name;}
    private:

};

class Exciton_Recombination : public Event{
    public:
        static const string name;
        string getName() const{return name;}
    private:
};

class Exciton_Dissociation : public Event{
    using Event::calculateExecutionTime;
    public:
        static const string name;
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,const double temperature,const double current_time){
            // Calculates dissociation using the Miller-Abrahams model
            double rate = prefactor*exp(-2*localization*distance);
            if(E_delta>0){
                rate *= exp(-E_delta/(K_b*temperature));
            }
            calculateExecutionTime(rate,current_time);
        }
        void calculateExecutionTime(const double prefactor,const double localization,const double distance,const double E_delta,const double reorganization,const double temperature,const double current_time){
            // Calculates dissociation using the Marcus model
            double rate = (prefactor/sqrt(4*Pi*reorganization*K_b*temperature))*exp(-2*localization*distance)*exp(-intpow(reorganization+E_delta,2)/(4*reorganization*K_b*temperature));
            calculateExecutionTime(rate,current_time);
        }
        string getName() const{return name;}
    private:
};

class Exciton_Intersystem_Crossing : public Event {
    public:
        static const string name;
        string getName() const{return name;}
    private:
};

class Exciton_Exciton_Annihilation : public Event{
    public:
        static const string name;
        string getName() const{return name;}
    private:
};

class Exciton_Polaron_Annihilation : public Event{
    public:
        static const string name;
        string getName() const{return name;}
    private:
};

#endif // EXCITON_H
