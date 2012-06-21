// All data use SI: dimension--m; density--Kg/m^3; pressure--Pa; time--second                  
#include "parameter.h"

namespace dem { 

///////////////////////////////////////////////////////////////////////////////////////
// Part A: These parameters do not change frequently
///////////////////////////////////////////////////////////////////////////////////////

// PI value
const REAL PI         = 3.14159265358979323846264338327950;

// gravitational acceleration
const REAL G          = 9.8;

// EPS (NOT float point relative precision, eps), problem domain dependent
const REAL EPS        = 1.0e-12;

// relative overlap between particles
const REAL MINOVERLAP = 1.0e-6;
const REAL MAXOVERLAP = 1.0e-2;

// measurable absolute overlap precision between particles, enabled/disabled by macro MEASURE_EPS
const REAL MEPS       = 1.0e-8;  // 0.1 micron or 0.01 micron

// random number seed
long idum             = -1;      // not a constant

// particle material property
const REAL YOUNG      = 2.90e+10;// quartz sand E  = 29GPa
const REAL POISSON    = 0.25;    // quartz sand v  = 0.25     
const REAL Gs         = 2.65;    // quartz sand Gs = 2.65    

// membrane particle material property
const REAL memYOUNG   = 1.40e+6; // 1.4MPa
const REAL memPOISSON = 0.49;

// other global variables
std::ofstream debugInf;          // print debugging information
int iteration;                   // iteration number

// output width and precision
const int OWID        = 16;      // 20, output width
const int OPREC       = 6;       // 10, output precision, number of digits after decimal dot

///////////////////////////////////////////////////////////////////////////////////////
// Part B: These parameters may change frequently and can be edited in main.cpp
///////////////////////////////////////////////////////////////////////////////////////

// number of MPI processes/partitions in each dimension
int NPX               = 1;
int NPY               = 1;
int NPZ               = 1;

// number of OpenMP threads
int  NUM_THREADS      = 1;       // default value, changed by command line

// 1. time integration method 
// dynamic
REAL TIMESTEP         = 5.0e-07; // time step
REAL MASS_SCL         = 1;       // mass scaling
REAL MNT_SCL          = 1;       // moment of inertia scaling
REAL GRVT_SCL         = 1;       // gravity scaling
REAL DMP_F            = 0;       // background viscous damping on mass   
REAL DMP_M            = 0;       // background viscous damping on moment of inertial
/*
// --- dynamic relaxation and scaling
dem::TIMESTEP         = 5.0e-06;
dem::MASS_SCL         = 1.0e+01;
dem::MNT_SCL          = 1.0e+01;
dem::GRVT_SCL         = 1.0e+03;
dem::DMP_F            = 2.0/dem::TIMESTEP;
dem::DMP_M            = 2.0/dem::TIMESTEP;
*/

// 2. normal damping and tangential friction
REAL DMP_CNT          = 0.30;    // damping ratio of viscous damping for normal contact force, for both particle-particle and particle-boundary contact
REAL FRICTION         = 0.50;    // constant coefficient of static friction between particles
REAL BDRYFRIC         = 0.50;    // constant coefficient of static friction between particle and rigid wall
REAL COHESION         = 0;       // 5.0e+8; cohesion between particles (10kPa)

// 3. boundary displacement rate
REAL COMPRESS_RATE    = 7.0e-03; // 7.0e-03 for triaxial; 1.0e-03 for isotropic and odometer.
REAL RELEASE_RATE     = 7.0e-03; // the same as above
REAL PILE_RATE        = 2.5e-01; // pile penetration velocity
REAL STRESS_ERROR     = 2.0e-02; // tolerance of stress equilibrium on rigid walls

} // namespace dem ends
