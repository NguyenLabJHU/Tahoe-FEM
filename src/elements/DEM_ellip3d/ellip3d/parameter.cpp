#include "parameter.h"

namespace dem { 

///////////////////////////////////////////////////////////////////////////////////////
// Part A: These parameters do not change frequently
///////////////////////////////////////////////////////////////////////////////////////

// PI value
const long double PI         = 3.141592653589793;

// gravitational acceleration
const long double G          = 9.8;

// numeric zero (NOT float point relative precision), problem domain dependent
const long double NUMZERO    = 1.0e-20;

// relative overlap between particles
const long double MINOVERLAP = 1.0e-6;
const long double MAXOVERLAP = 1.0e-2;

// random number seed
long idum                    = -1;      // not a constant

// particle material property
const long double YOUNG      = 2.90e+10;// quartz sanxd E=29GPa
const long double POISSON    = 0.25;    // quartz sand v= 0.25     
const long double Gs         = 2.65;    // quartz sand Gs= 2.65    

// other global variables
std::ofstream g_debuginf;               // print debugging information
int g_iteration;                        // iteration number

// output width and precision
const int OWID               = 20;      // output width
const int OPREC              = 10;      // output precision, number of digits after decimal dot 

///////////////////////////////////////////////////////////////////////////////////////
// Part B: These parameters may change frequently and can be easily edited in main.cpp
///////////////////////////////////////////////////////////////////////////////////////

// number of OpenMP threads
int         NUM_THREADS      = 1;

// 1. time integration method 
long double TIMESTEP         = 5.0e-07; // time step
long double MASS_SCL         = 1;       // mass scaling
long double MNT_SCL          = 1;       // moment of inertia scaling
long double GRVT_SCL         = 1;       // gravity scaling
long double DMP_F            = 0;       // background viscous damping on mass   
long double DMP_M            = 0;       // background viscous damping on moment of inertial

// 2. normal damping and tangential friction
long double DMP_CNT          = 0.05;    // damping ratio of viscous damping for normal contact force, for both particle-particle and particle-boundary contact
long double FRICTION         = 0.5;     // constant coefficient of static friction between particles
long double BDRYFRIC         = 0.5;     // constant coefficient of static friction between particle and rigid wall
long double COHESION         = 5.0e+8;  // cohesion between particles (10kPa)

// 3. boundary displacement rate
long double COMPRESS_RATE    = 7.0e-03; // 7.0e-03 for triaxial; 1.0e-03 for isotropic and odometer.
long double RELEASE_RATE     = 7.0e-03; // the same as above
long double PILE_RATE        = 2.5e-01; // pile penetration velocity
long double STRESS_ERROR     = 2.0e-02; // tolerance of stress equilibrium on rigid walls

} // namespace dem ends
