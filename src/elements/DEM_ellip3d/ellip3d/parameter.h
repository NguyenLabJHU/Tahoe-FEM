// These paramenters may need to change frequently and it should be done in main.cpp

// All data use SI: dimension--m; density--Kg/m^3; pressure--Pa; time--second                  

#ifndef PARAMETER_H
#define PARAMETER_H

namespace dem { 

////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1. time integration method 
long double TIMESTEP      = 5.0e-07; // time step
long double MASS_SCL      = 1;       // mass scaling
long double MNT_SCL       = 1;       // moment of inertia scaling
long double GRVT_SCL      = 1;       // gravity scaling
long double DMP_F         = 0;       // background viscous damping on mass   
long double DMP_M         = 0;       // background viscous damping on moment of inertial
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2. normal damping and tangential friction
long double DMP_CNT       = 0.05;    // damping ratio of viscous damping for normal contact force, for both particle-particle and particle-boundary contact
long double FRICTION      = 0.5;     // constant coefficient of static friction between particles
long double BDRYFRIC      = 0.5;     // constant coefficient of static friction between particle and rigid wall
long double COHESION      = 5.0e+8;  // cohesion between particles (10kPa)
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3. boundary displacement rate
long double COMPRESS_RATE = 7.0e-03; // 7.0e-03 for triaxial; 1.0e-03 for isotropic and odometer.
long double RELEASE_RATE  = 7.0e-03;
long double PILE_RATE     = 2.5e-01; // pile penetration velocity
long double STRESS_ERROR  = 2.0e-02; // tolerance of stress equilibrium on rigid walls
////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// 4. other global variables
long     idum             = -1;      // random number seed
std::ofstream g_exceptioninf;        // record debugging information
int      g_iteration;                // iteration number 
////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace dem ends

#endif
