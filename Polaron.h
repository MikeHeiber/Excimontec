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

using namespace std;

class Polaron : public Object{
    public:
        static const string name;
        Polaron(const double time,const int tag_num,const Coords& start_coords) : Object(time,tag_num,start_coords){}
        string getName(){return name;}
    private:

};

class Polaron_Hop : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature, const double prefactor){
            setDestCoords(dest_coords);
            // No target object
            // setWaitTime
        }
        string getName(){return name;}
    private:

};

class Electron_Hole_Recombination : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature, const double prefactor){
            setDestCoords(dest_coords);
            // No target object
            // setWaitTime
        }
        string getName(){return name;}
    private:

};

class Polaron_Extraction : public Event{
    public:
        static const string name;
        void calculateEvent(const Coords& dest_coords,const double distance,const double E_delta,const int temperature, const double prefactor){
            setDestCoords(dest_coords);
            // No target object
            // setWaitTime
        }
        string getName(){return name;}
    private:
};

#endif // POLARON_H
