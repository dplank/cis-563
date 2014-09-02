// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "vec.h"

#define LERP(a,b,t) (1-t)*a + t*b


// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;
extern const double S;
extern const double rho_f;
extern const double rho_h;

extern const double Tmax;     
extern const double Tignition;  
extern const double CoolT;    
extern const double Radius;
extern const vec3 center;

#endif