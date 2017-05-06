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
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            // No destination coords.  Destination coords are chosen upon execution.
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:


};

class Exciton_Hop : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            setDestCoords(dest_coords);
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:

};

class Exciton_Recombination : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            // No destination site
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:
};

class Exciton_Dissociation : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            setDestCoords(dest_coords);
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:
};

class Exciton_Intersystem_Crossing : public Event {
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            // No destination coords
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:
};

class Exciton_Exciton_Annihilation : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            setDestCoords(dest_coords);
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:
};

class Exciton_Polaron_Annihilation : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double rate,const double current_time){
            setDestCoords(dest_coords);
            // No target object
            setExecutionTime(current_time+(-1/rate)*log(rand01()));
        }
        string getName() const{return name;}
    private:
};

#endif // EXCITON_H
