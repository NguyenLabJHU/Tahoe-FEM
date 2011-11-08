#ifndef PARAMETER_H
#define PARAMETER_H

#include <fstream>

namespace dem { 

///////////////////////////////////////////////////////////////////////////////////////
// Part A: These parameters do not change frequently
///////////////////////////////////////////////////////////////////////////////////////

// PI value
extern const long double PI;

// Gravitational acceleration
extern const long double G;

// absolute numeric precision
extern const long double PREC;

// random number seed
extern long idum;

// random shape for each particle
//#define RANDOM_SHAPE

// particle material property
extern const long double YOUNG;  
extern const long double POISSON;      
extern const long double Gs;     

// other global variables
extern std::ofstream g_exceptioninf;
extern int g_iteration;

///////////////////////////////////////////////////////////////////////////////////////
// Part B: These parameters may change frequently and can be easily edited in main.cpp
///////////////////////////////////////////////////////////////////////////////////////

// number of OpenMP threads
extern int         NUM_THREADS;

// 1. time integration method 
extern long double TIMESTEP;
extern long double MASS_SCL;
extern long double MNT_SCL;
extern long double GRVT_SCL;
extern long double DMP_F;
extern long double DMP_M;

// 2. normal damping and tangential friction
extern long double DMP_CNT;
extern long double FRICTION;
extern long double BDRYFRIC;
extern long double COHESION;

// 3. boundary displacement rate
extern long double COMPRESS_RATE;
extern long double RELEASE_RATE;
extern long double PILE_RATE;
extern long double STRESS_ERROR;

} // namespace dem ends

#endif
