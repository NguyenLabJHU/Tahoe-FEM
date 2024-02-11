// #define marcos are placed in namespace dem nominally for convenience,
// they are not contained by any namespace, i.e., they are preprocessed.
 
#ifndef PARAMETER_H
#define PARAMETER_H
#include "realtypes.h"
#include <fstream>

namespace dem { 

///////////////////////////////////////////////////////////////////////////////////////
// Part A: These parameters do not change frequently
///////////////////////////////////////////////////////////////////////////////////////

// PI value
extern const REAL PI;

// Gravitational acceleration
extern const REAL G;

// EPS (NOT float point relative precision, eps), problem domain dependent
extern const REAL EPS;

// relative overlap between particles
extern const REAL MINOVERLAP;
extern const REAL MAXOVERLAP;

// macro to toggle on/off MEPS
#define MEASURE_EPS

// measurable absolute overlap precision between particles
extern const REAL MEPS;

// random number seed
extern long idum;

// macro to toggle on/off random shape for each particle
//#define RANDOM_SHAPE

// particle material property
extern const REAL YOUNG;  
extern const REAL POISSON;      
extern const REAL Gs;   

// critical tensile stress for particle sub-division
extern const REAL sigmaCritical;	// calculated from experiment

// compressive strength for particle sub-division based on Hoek-Brown criterion
extern const REAL sigmaCompress;	// calculated from experiment
extern const REAL mi;			// material const, for granite mi=32.4

// critical maximum tensile stress for contact point criterion
extern const REAL ContactTensileCritical;	// calculated from experiment

// Weibull modulus used for particle strength
extern const REAL weibullModulus;
extern const REAL basicRadius;	// the radius of the base particle in weibull function

// properties for the springs
extern const REAL sigma_f;	// soft criterion for spring, need to be calibrated
extern const REAL Cf;		// crack propagate speed, calculated from experiment

// membrane particle material property
extern const REAL memYOUNG;
extern const REAL memPOISSON;

// other global variables
extern std::ofstream g_debuginf;
extern std::ofstream g_timeinf;
extern int g_iteration;

// number of particles broken in different subdivisions
extern int numBrokenType1;	// Hoek-Brown criterion or maximum shear stress
extern int numBrokenType2;	// maximum tensile stress at contacts

// output field width and precision
extern const int OWID;
extern const int OPREC;

// number of timesteps for transition process
extern const int numStepTransition;

///////////////////////////////////////////////////////////////////////////////////////
// Part B: These parameters may change frequently and can be easily edited in main.cpp
///////////////////////////////////////////////////////////////////////////////////////

// number of OpenMP threads
extern int  NUM_THREADS;

// 1. time integration method 
extern REAL TIMESTEP;
extern REAL MASS_SCL;
extern REAL MNT_SCL;
extern REAL GRVT_SCL;
extern REAL DMP_F;
extern REAL DMP_M;

// 2. normal damping and tangential friction
extern REAL DMP_CNT;
extern REAL FRICTION;
extern REAL BDRYFRIC;
extern REAL COHESION;

// 3. boundary displacement rate
extern REAL COMPRESS_RATE;
extern REAL RELEASE_RATE;
extern REAL PILE_RATE;
extern REAL STRESS_ERROR;

} // namespace dem ends

#endif
