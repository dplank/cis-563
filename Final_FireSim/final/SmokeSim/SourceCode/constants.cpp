// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {20, 20, 9};
#endif

const double theCellSize = 0.5;
const double S = 0.1;
const double rho_f = 1.0;
const double rho_h= 0.1;

const double Tmax = 50.0;   
const double Tignition = 20.0; 
const double CoolT = 2.0;  

const double Radius = 2.5 * theCellSize + 0.1;
const vec3 center = vec3(10,5,5);