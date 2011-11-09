////////////////////////////////////////////////////////////////////////////////////////////////////
// File: main.cpp, the main program, controlling simulation type and parameters  
// All data use SI: dimension--m; density--Kg/m^3; pressure--Pa; time--second                  
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// header files
#include "parameter.h"
#include "gradation.h"
#include "rectangle.h"
#include "assembly.h"
#include <time.h>
#include <iostream>
#include <vector>
#include <cstdlib>

////////////////////////////////////////////////////////////////////////////////////////////////////
// main Program
int main(int argc, char* argv[])
{
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Part 0: command line arguments and timestamps
    time_t time1, time2;
    time(&time1);
    if (argc == 2) dem::NUM_THREADS = atoi(argv[1]);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Part 1: setup parameters (override parameter.cpp)
    // 1. time integration method
    // --- dynamic
    /*
    dem::TIMESTEP      = 5.0e-07; // time step
    dem::MASS_SCL      = 1;       // mass scaling
    dem::MNT_SCL       = 1;       // moment of inertial scaling
    dem::GRVT_SCL      = 1;       // gravity scaling
    dem::DMP_F         = 0;       // background viscous damping on mass
    dem::DMP_M         = 0;       // background viscous damping on moment of inertial
    */

    // --- dynamic relaxation and scaling
    dem::TIMESTEP      = 5.0e-06;
    dem::MASS_SCL      = 1.0e+01;
    dem::MNT_SCL       = 1.0e+01;
    dem::GRVT_SCL      = 0;       // 1.0e+03;
    dem::DMP_F         = 2.0/dem::TIMESTEP;
    dem::DMP_M         = 2.0/dem::TIMESTEP;

    // 2. normal damping and tangential friciton
    dem::DMP_CNT       = 0.05;    // normal contact damping ratio
    dem::FRICTION      = 0.5;     // coefficient of friction between particles
    dem::BDRYFRIC      = 0.5;     // coefficient of friction between particle and rigid wall
    dem::COHESION      = 0;//5.0e+8;  // cohesion between particles

    // 3. boundary displacement rate
    dem::COMPRESS_RATE = 7.0e-03; // 7.0e-03 for triaxial; 1.0e-03 for isotropic and odometer.
    dem::RELEASE_RATE  = 7.0e-03; // the same as above
    dem::PILE_RATE     = 2.5e-01; // pile penetration velocity
    dem::STRESS_ERROR  = 2.0e-02; // tolerance of stress equilibrium on rigid walls

    dem::assembly A;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Part 2: set up a simulation to run
    A.triaxial(100000,             // total_steps
	       100,                // number of snapshots
	       10,                 // record interval
	       5.0e+5,             // pre-existing ambient pressure from initial isotropic compression
	       "iso_particle_500k",// input file, initial particles
	       "iso_boundary_500k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary",     // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Part 3: record run time
    time(&time2);
    std::cout << std::endl 
	      << "simulation start time: "<<ctime(&time1);
    std::cout << "simulation  end  time: "<<ctime(&time2);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
// Notes:
// freetype (setting of free particles):
// 0 - a single free particle
// 1 - a layer of free particles
// 2 - multiple layers of free particles

// particle type (setting for each particle):
//  0 - free particle
//  1 - fixed particle
//  2 - special case 2 (pure moment): translate first, then rotate only, MNT_START needs to be defined
//  3 - special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
//  4 - special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only
//  5 - free boundary particle
// 10 - ghost particle

////////////////////////////////////////////////////////////////////////////////////////////////////
// Various types of simulation, copy them into part 2 of main() to run, only one block a time.
/*

    // for triaxalPtclBdry
    A.TrimPtclBdryByHeight(0.061,
			   "dep_particle_end",
			   "dep_particle_trimmed");

    // for de-fe coupling			   
    A.TrimPtclBdryByHeight(0.061,
			   "pile_particle_ini",
			   "pile_particle_trimmed");

    A.triaxialPtclBdryIni(100000,             // total_steps
			  100,                // number of snapshots
                          10,                 // record interval
			  2.5e+6,             // confining pressure to achieve
			  "ini_particle_500k",// input file, initial particles
			  "ini_boundary_500k",// input file, initial boundaries
			  "ini_particle",     // output file, resulted particles, including snapshots 
			  "ini_boundary",     // output file, resulted boundaries
			  "ini_contact",      // output file, resulted contacts, including snapshots 
			  "ini_progress",     // output file, statistical info
			  "ini_debug");       // output file, debug info

    A.triaxialPtclBdryIni(100000,             // total_steps
			  100,                // number of snapshots
                          10,                 // record interval
			  5.0e+5,             // confining pressure to achieve
			  "ini_particle_ini", // input file, initial particles
			  "ini_boundary_ini", // input file, initial boundaries
			  "ini_particle",     // output file, resulted particles, including snapshots 
			  "ini_boundary",     // output file, resulted boundaries
			  "ini_contact",      // output file, resulted contacts, including snapshots 
			  "ini_progress",     // output file, statistical info
			  "ini_debug");       // output file, debug info


    A.triaxialPtclBdry(100000,             // total_steps
		       100,                // number of snapshots
                       10,                 // record interval
		       "tri_particle_ini", // input file, initial particles
		       "tri_boundary_ini", // input file, initial boundaries
		       "tri_particle",     // output file, resulted particles, including snapshots 
		       "tri_boundary",     // output file, resulted boundaries
		       "tri_contact",      // output file, resulted contacts, including snapshots 
		       "tri_progress",     // output file, statistical info
		       "tri_balanced",     // output file, balanced status
		       "tri_debug");       // output file, debug info


    // size, shape, and gradation of particles
    int rorc             = 1;     // rectangular = 1 or cylindrical = 0
    long double dimn     = 0.05;  // specimen dimension
    long double ratio_ba = 0.8;   // ratio of radius b to radius a
    long double ratio_ca = 0.6;   // ratio of radius c to radius a
    std::vector<long double> percent;  // mass percentage of particles smaller than a certain size
    std::vector<long double> ptclsize; // particle size
    percent.push_back(1.00); ptclsize.push_back(2.5e-3);
    //percent.push_back(0.80); ptclsize.push_back(2.3e-3);
    //percent.push_back(0.60); ptclsize.push_back(2.0e-3);
    //percent.push_back(0.30); ptclsize.push_back(1.5e-3);
    //percent.push_back(0.10); ptclsize.push_back(1.0e-3);
    dem::gradation grad(rorc, dimn, ratio_ba, ratio_ca, percent.size(), percent, ptclsize);
    A.deposit_RgdBdry(grad,
		      2,                  // freetype, setting of free particles 
		      100000,             // total_steps
		      100,                // number of snapshots
                      10,                 // record interval
		      3.0,                // height of floating particles relative to dimn
		      "flo_particle_end", // output file, initial particles
		      "dep_boundary_ini", // output file, initial boundaries
		      "dep_particle",     // output file, resulted particles, including snapshots 
		      "dep_contact",      // output file, resulted contacts, including snapshots 
		      "dep_progress",     // output file, statistical info
		      "cre_particle",     // output file, resulted particles after trmming
		      "cre_boundary",     // output file, resulted boundaries after trmming
		      "dep_debug");       // output file, debug info


    // size, shape, and gradation of particles
    int rorc             = 1;     // rectangular = 1 or cylindrical = 0
    long double dimn     = 0.05;  // specimen dimension
    long double ratio_ba = 0.8;   // ratio of radius b to radius a
    long double ratio_ca = 0.6;   // ratio of radius c to radius a
    std::vector<long double> percent;  // mass percentage of particles smaller than a certain size
    std::vector<long double> ptclsize; // particle size
    percent.push_back(1.00); ptclsize.push_back(2.5e-3);
    //percent.push_back(0.80); ptclsize.push_back(2.0e-3);
    //percent.push_back(0.60); ptclsize.push_back(1.6e-3);
    //percent.push_back(0.30); ptclsize.push_back(1.0e-3);
    //percent.push_back(0.10); ptclsize.push_back(0.5e-3);
    dem::gradation grad(rorc, dimn, ratio_ba, ratio_ca, percent.size(), percent, ptclsize);
    A.deposit_PtclBdry(grad,
		       2,                  // freetype, setting of free particles
		       1.0,                // relative container size, 0.8/1.0/1.2---small/medium/large
		       100000,             // total_steps
		       100,                // number of snapshots
                       10,                 // record interval
		       "flo_particle_end", // output file, initial particles
		       "dep_particle",     // output file, resulted particles, including snapshots 
		       "dep_contact",      // output file, resulted contacts, including snapshots 
		       "dep_progress",     // output file, statistical info
		       "dep_debug");       // output file, debug info
 
   
    A.scale_PtclBdry(20000,             // total_steps
		     100,               // number of snapshots  
                     10,                // record interval
		     0.05,              // dimension of particle-composed-boundary
		     1.0,               // relative container size, 0.8/1.0/1.2---small/medium/large
		     "dep_particle_end",// input file, initial particles
		     "scl_particle",    // output file, resulted particles, including snapshots 
		     "scl_contact",     // output file, resulted contacts, including snapshots
		     "scl_progress",    // output file, statistical info
		     "scl_debug");      // output file, debug info


    A.ellipPile_Disp(50000,              // total_steps
		     100,                // number of snapshots
                     10,                 // record interval
		     0.05,               // dimension of particle-composed-boundary
		     1.0,                // relative container size, 0.8/1.0/1.2---small/medium/large
		     "pile_particle_ini",// input file, initial particles, an ellipsoidal pile info added
		     "pile_particle",    // output file, resulted particles, including snapshots 
		     "pile_contact",     // output file, resulted contacts, including snapshots 
		     "pile_progress",    // output file, statistical info
		     "pile_debug");      // output file, debug info


    A.rectPile_Disp(50000,              // total_steps
		    100,                // number of snapshots
                    10,                 // record interval
		    "pile_particle_ini",// input file, initial particles
		    "pile_boundary_ini",// input file, initial boundaries, rectangular pile boundary info added
		    "pile_particle",    // output file, resulted particles, including snapshots 
		    "pile_boundary",    // output file, resulted boundaries
		    "pile_contact",     // output file, resulted contacts, including snapshots 
		    "pile_progress",    // output file, statistical info
		    "pile_debug");      // output file, debug info


    A.ellipPile_Impact(50000,              // total_steps
		       100,                // number of snapshots
                       10,                 // record interval
		       0.05,               // size of particle-composed-boundary
		       "ipt_particle_ini", // input file, initial particles, an ellipsoidal pile info added
		       "dep_boundary_ini", // input file, initial boundaries
		       "ipt_particle",     // output file, resulted particles, including snapshots 
		       "ipt_contact",      // output file, resulted contacts, including snapshots 
		       "ipt_progress",     // output file, statistical info
		       "ipt_debug");       // output file, debug info


    A.ellipPile_Impact_p(50000,              // total_steps
			 100,                // number of snapshots
                         10,                 // record interval
			 0.05,               // size of particle-composed-boundary
			 "ipt_particle_ini", // input file, initial particles, an ellipsoidal pile info added
			 "ipt_particle",     // output file, resulted particles, including snapshots 
			 "ipt_contact",      // output file, resulted contacts, including snapshots 
			 "ipt_progress",     // output file, statistical info
			 "ipt_debug");       // output file, debug info

    
    A.deposit(5000,               // total_steps
	      100,                // number of snapshots
              10,                 // record interval
	      "flo_particle_end", // input file, initial particles
	      "dep_boundary_ini", // input file, initial boundaries
	      "dep_particle",     // output file, resulted particles, including snapshots 
	      "dep_contact",      // output file, resulted contacts, including snapshots 
	      "dep_progress",     // output file, statistical info
	      "dep_debug");       // output file, debug info


    A.squeeze(300000,             // total_steps
	      100000,             // initial_steps to reach equilibrium
              100,                // number of snapshots
              10,                 // record interval
	      +1,                 // -1 squeeze; +1 loosen
	      "flo_particle_end", // input file, initial particles
	      "dep_boundary_ini", // input file, initial boundaries
	      "dep_particle",     // output file, resulted particles, including snapshots 
	      "dep_boundary",     // output file, resulted boundaries
	      "dep_contact",      // output file, resulted contacts, including snapshots 
	      "dep_progress",     // output file, statistical info
	      "dep_debug");       // output file, debug info


    A.collapse(rorc,
	       100000,
	       100,
               10,                // record interval
	       "cre_particle",    // input file, initial particles
	       "clp_boundary",    // output file, initial boundaries
	       "clp_particle",    // output file, resulted particles, including snapshots
	       "clp_contact",     // output file, resulted contacts, including snapshots 
	       "clp_progress",    // output file, statistical info
	       "clp_debug");      // output file, debug info


    A.isotropic(100000,             // total_steps
		100,                // number of snapshots
                10,                 // record interval
		1.0e+3,             // a low confining pressure to achieve for initialization
		"cre_particle",     // input file, initial particles
		"cre_boundary",     // input file, initial boundaries
		"iso_particle",     // output file, resulted particles, including snapshots 
		"iso_boundary",     // output file, resulted boundaries 
		"iso_contact",      // output file, resulted contacts, including snapshots 
		"iso_progress",     // output file, statistical info
		"iso_balanced",     // output file, balanced status
		"iso_debug");       // output file, debug info

    
    A.isotropic(100000,             // total_steps
		100,                // number of snapshots
		1.0e+3,             // pre-existing confining pressure from initial isotropic compression
		1.0e+5,             // confining pressure for final isotropic compression
		100,                // fine steps for applying pressure
		"iso_particle_1k",  // input file, initial particles
		"iso_boundary_1k",  // input file, initial boundaries
		"iso_particle",     // output file, resulted particles, including snapshots 
		"iso_boundary",     // output file, resulted boundaries
		"iso_contact",      // output file, resulted contacts, including snapshots 
		"iso_progress",     // output file, statistical info
		"iso_balanced",     // output file, balanced status
		"iso_debug");       // output file, debug info

    
    A.isotropic(100000,             // total_steps
		100,                // number of snapshots
                10,                 // record interval
		1.0e+5,             // pre-existing confining pressure from initial isotropic compression
		1.5e+6,             // confining pressure for final isotropic compression
		100,                // fine steps for applying pressure
		"iso_particle_100k",// input file, initial particles
		"iso_boundary_100k",// input file, initial boundaries
		"iso_particle",     // output file, resulted particles, including snapshots 
		"iso_boundary",     // output file, resulted boundaries
		"iso_contact",      // output file, resulted contacts, including snapshots 
		"iso_progress",     // output file, statistical info
		"iso_balanced",     // output file, balanced status
		"iso_debug");       // output file, debug info

    
    long double sigma_values[4]={1.0e+5, 5.0e+5, 1.0e+5, 7.0e+5}; // last one must be a larger value	
    A.isotropic(100000,             // total_steps
		100,                // number of snapshots
                10,                 // record interval
		4,                  // number of points indicating pressure applying process
		sigma_values,       // loading, unloading and reloading stress path
		100,                // fine steps for applying pressure between two adjacent values
		"iso_particle_100k",// input file, initial particles
		"iso_boundary_100k",// input file, initial boundaries
		"iso_particle",     // output file, resulted particles, including snapshots 
		"iso_boundary",     // output file, resulted boundaries
		"iso_contact",      // output file, resulted contacts, including snapshots 
		"iso_progress",     // output file, statistical info
		"iso_balanced",     // output file, balanced status
		"iso_debug");       // output file, debug info

    
    A.odometer(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // record interval
	       1.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       1.5e+6,             // major principle stress for final odometer compression
	       100,                // fine steps for applying pressure
	       "iso_particle_100k",// input file, initial particles
	       "iso_boundary_100k",// input file, initial boundaries
	       "odo_particle",     // output file, resulted particles, including snapshots 
	       "odo_boundary",     // output file, resulted boundaries
	       "odo_contact",      // output file, resulted contacts, including snapshots 
	       "odo_progress",     // output file, statistical info
	       "odo_balanced",     // output file, progress odometer balanced status
	       "odo_debug");       // output file, debug info

    
    long double sigma_values[4]={1.0e+5, 5.0e+5, 1.0e+5, 1.0e+6}; // last one must be a larger value	
    A.odometer(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // record interval
	       4,                  // number of points indicating pressure applying process
	       sigma_values,       // loading, unloading and reloading stress path
	       100,                // fine steps for applying pressure
	       "iso_particle_100k",// input file, initial particles
	       "iso_boundary_100k",// input file, initial boundaries
	       "odo_particle",     // output file, resulted particles, including snapshots 
	       "odo_boundary",     // output file, resulted boundaries
	       "odo_contact",      // output file, resulted contacts, including snapshots 
	       "odo_progress",     // output file, statistical info
	       "odo_balanced",     // output file, progress odometer balanced status
	       "odo_debug");       // output file, debug info


    A.triaxial(100000,             // total_steps
	       100,                // number of snapshots
	       1.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_100k",// input file, initial particles
	       "iso_boundary_100k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary",     // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    A.triaxial(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // record interval
	       3.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_300k",// input file, initial particles
	       "iso_boundary_300k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary"      // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    A.triaxial(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // record interval
	       5.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_500k",// input file, initial particles
	       "iso_boundary_500k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary",     // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    A.triaxial(120000,             // total_steps
	       30000,              // at which step to unload
	       100,                // number of snapshots
               10,                 // record interval
	       3.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_300k",// input file, initial particles
	       "iso_boundary_300k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary",     // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    A.truetriaxial(100000,             // total_steps
		   100,                // number of snapshots
                   10,                 // record interval
		   1.0e+4,             // pre-existing confining pressure from initial isotropic compression
		   1.0e+5,             // sigma_w
		   3.0e+5,             // sigma_l
		   9.0e+5,             // sigma_h
		   100,                // fine steps for applying pressure
		   "iso_particle_10k", // input file, initial particles
		   "iso_boundary_10k", // input file, initial boundaries
		   "tru_particle",     // output file, resulted particles, including snapshots 
		   "tru_boundary",     // output file, resulted boundaries
		   "tru_contact",      // output file, resulted contacts, including snapshots 
		   "tru_progress",     // output file, statistical info
		   "tru_balanced",     // output file, balanced status
		   "tru_debug");       // output file, debug info


    A.unconfined(100000,             // total_steps
		 100,                // number of snapshots
                 10,                 // record interval
		 "flo_particle_end", // input file, initial particles
		 "unc_boundary",     // input file, initial boundaries
		 "unc_particle",     // output file, resulted particles, including snapshots 
		 "unc_contact",      // output file, resulted contacts, including snapshots 
		 "unc_progress",     // output file, statistical info
		 "unc_debug");       // output file, debug info
*/
