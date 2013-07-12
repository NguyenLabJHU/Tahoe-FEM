#include "realtypes.h"
#include "const.h"
#include "Parameter.h"
#include "Gradation.h"
#include "Rectangle.h"
#include "Assembly.h"
#include "Vec.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/timer/timer.hpp>

// serialization of pointers to objects of derived classes
// http://www.boost.org/doc/libs/1_50_0/libs/serialization/doc/serialization.html#registration
// this is called registration, and this method of registration is referred to as "key export".
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(dem::plnBoundary, "plnBoundary")

// optimization on non-template types which have a fixed amount of data stored at fixed field positions
BOOST_IS_MPI_DATATYPE(dem::Vec)
BOOST_IS_MPI_DATATYPE(dem::Particle)
BOOST_IS_MPI_DATATYPE(dem::Contact)

int main(int argc, char* argv[]) {

  boost::mpi::environment  boostEnv(argc, argv);
  boost::mpi::communicator boostWorld;

  boost::timer::auto_cpu_timer boostTimer;

  if (boostWorld.rank() == 0) {
    if (argc != 2) {
      std::cout << "please specify data file in the form: paraEllip3d input.txt" << std::endl;
      return -1;  
    }

    dem::Parameter::getSingleton().readIn(argv[1]);
    //dem::Parameter::getSingleton().writeOut();
    int mpiProcX = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcX"]);
    int mpiProcY = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcY"]);
    int mpiProcZ = static_cast<int> (dem::Parameter::getSingleton().parameter["mpiProcZ"]);
    if (mpiProcX * mpiProcY * mpiProcZ != boostWorld.size() ) {
      std::cout << "number of MPI processes does match grids in data file!" << std::endl;
      return -1;
    }
  }

  broadcast(boostWorld, dem::Parameter::getSingleton(), 0); // broadcast from root process 0

  dem::debugInf.open("debugInf");
  if(!dem::debugInf) { std::cout << "stream error!" << std::endl; exit(-1);}
  dem::debugInf.setf(std::ios::scientific, std::ios::floatfield);

  dem::Assembly assemb;
  assemb.setCommunicator(boostWorld);

  int simuType = static_cast<int> (dem::Parameter::getSingleton().parameter["simuType"]);
  switch (simuType) {
  case 1: // deposit spatially scattered particles into a rigid container
    assemb.depositIntoContainer();
    break;
  case 2: // resume deposition using specified data file of particles and boundaries
    assemb.resumeDepositIntoContainer(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
				      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    break;
  case 3: // expand particles inside a virtual cavity and see what occurs
    assemb.expandCavityParticle();
    break;
  case 4: // resume expanding particles inside a virtual cavity and see what occurs
    assemb.resumeExpandCavityParticle(dem::Parameter::getSingleton().datafile["boundaryFile"].c_str(),
				      dem::Parameter::getSingleton().datafile["particleFile"].c_str());
    break;   
  }
  
  dem::debugInf.close();

  ///////////////////////////////////////////////////////////////////////////////
  // Part 2: set up a simulation to run

  /*
  // container properties
  dem::Rectangle container(0.05, 0.05, 0.05, vec(0, 0, 0)); // dimx, dimy, dimz, center
  // particle shape, size and percentage
  REAL ptcl_ratio_ba = 0.8;  // ratio of radius b to radius a
  REAL ptcl_ratio_ca = 0.6;  // ratio of radius c to radius a
  std::vector<REAL> percent; // mass percentage of particles smaller than a certain size
  std::vector<REAL> ptclSize;    // particle size
  percent.push_back(1.00); ptclSize.push_back(2.5e-3);
  percent.push_back(0.80); ptclSize.push_back(2.3e-3);
  percent.push_back(0.60); ptclSize.push_back(2.0e-3);
  percent.push_back(0.30); ptclSize.push_back(1.5e-3);
  percent.push_back(0.10); ptclSize.push_back(1.0e-3);
  dem::Gradation ptclGradation(percent.size(), percent, ptclSize, ptcl_ratio_ba, ptcl_ratio_ca);
  assemb.setContainer(container);
  assemb.setGradation(ptclGradation);
  */

  /*
  assemb.plotBoundary("container.plot");
  assemb.depositIntoContainer
    (2,                  // particleLayers, setting of free particles 
     1000000,            // total_steps
     100,                // number of snapshots
     10,                 // print interval
     5.0,                // relative height of floating particles based on container height
     "flo_particle_end", // output file, initial particles for depositing
     "dep_boundary_ini", // output file, initial boundaries for depositing
     "dep_particle",     // output file, resulted particles, including snapshots 
     "dep_contact",      // output file, resulted contacts, including snapshots 
     "dep_progress",     // output file, statistical info
     "trm_particle_end", // output file, resulted particles after trmming
     "trm_boundary_end", // output file, resulted boundaries after trmming
     "dep_debug");       // output file, debug 
  */

  /*
  assemb.deposit(100,                // total_steps
		 10,                  // number of snapshots
		 10,                  // print interval
		 "dep_particle_end", // input file, initial particles
		 "dep_boundary_ini", // input file, initial boundaries
		 "exp_particle",     // output file, resulted particles, including snapshots 
		 "exp_contact",      // output file, resulted contacts, including snapshots 
		 "exp_progress",     // output file, statistical info
		 "exp_debug");       // output file, debug info
  */
  /*
  dem::Rectangle cavity(0.05*dimx, 0.05*dimy, 0.05*dimz, dem::vec(0, 0, -0.25*dimz));
  assemb.setCavity(cavity);
  assemb.plotBoundary("container.plot");
  assemb.plotCavity("cavity.plot");
  assemb.expandCavityParticles(true, 
			       0.40,
			       "cav_particle_ini",
			       "dep_particle_end", 
			       "exp_particle_ini");
  assemb.deposit(1000000,            // total_steps
		 100,                // number of snapshots
		 10,                 // print interval
		 "exp_particle_ini", // input file, initial particles
		 "dep_boundary_ini", // input file, initial boundaries
		 "exp_particle",     // output file, resulted particles, including snapshots 
		 "exp_contact",      // output file, resulted contacts, including snapshots 
		 "exp_progress",     // output file, statistical info
		 "exp_debug");       // output file, debug info  
  */
  /*
  dem::Rectangle cavity(0.1*dimx, 0.1*dimy, 0.1*dimz, dem::vec(0, 0, -0.25*dimz));
  assemb.setCavity(cavity);
  assemb.trimCavity(true, 
		    "dep_particle_end",
		    "cav_particle_ini");
  assemb.buildCavityBoundary(6, "cav_boundary_ini"); // existMaxId; output file
  assemb.deposit(250000,             // total_steps
		 100,                // number of snapshots
		 10,                 // print interval
		 "cav_particle_ini", // input file, initial particles
		 "dep_boundary_ini", // input file, initial boundaries
		 "cav_particle",     // output file, resulted particles, including snapshots 
		 "cav_contact",      // output file, resulted contacts, including snapshots 
		 "cav_progress",     // output file, statistical info
		 "cav_debug");       // output file, debug info 
  */
  /*
  dem::Rectangle cavity(0.1*dimx, 0.1*dimy, 0.1*dimz, dem::vec(0, 0, -0.25*dimz));
  assemb.setCavity(cavity);
  assemb.trimCavity(true, 
		    "dep_particle_end",
		    "cav_particle_ini");
  assemb.buildCavityBoundary(6, "cav_boundary_ini"); // existMaxId; output file
  assemb.plotBoundary("container.plot");
  assemb.plotCavity("cavity.plot");
  assemb.depositAfterCavity(250000,             // total_steps
			    100,                // number of snapshots
			    10,                 // print interval
			    "cav_particle_ini", // input file, initial particles
			    "dep_boundary_ini", // input file, initial boundaries
			    "cav_boundary_ini", // input file, initial cavity boundaries
			    "cav_particle",     // output file, resulted particles, including snapshots 
			    "cav_contact",      // output file, resulted contacts, including snapshots 
			    "cav_progress",     // output file, statistical info
			    "cav_debug");       // output file, debug info 
  */
  /*
  assemb.deposit(1,                  // total_steps
		 1,                  // number of snapshots
		 1,                  // print interval
		 "dep_particle_end", // input file, initial particles
		 "dep_boundary_ini", // input file, initial boundaries
		 "tst_particle",     // output file, resulted particles, including snapshots 
		 "tst_contact",      // output file, resulted contacts, including snapshots 
		 "tst_progress",     // output file, statistical info
		 "tst_debug");       // output file, debug info  
  */
  /*
  assemb.depositIntoContainer
    (2,                  // particleLayers, setting of free particles 
     1000000,            // total_steps
     100,                // number of snapshots
     10,                 // print interval
     5.0,                // relative height of floating particles based on container height
     "flo_particle_end", // output file, initial particles for depositing
     "dep_boundary_ini", // output file, initial boundaries for depositing
     "dep_particle",     // output file, resulted particles, including snapshots 
     "dep_contact",      // output file, resulted contacts, including snapshots 
     "dep_progress",     // output file, statistical info
     "trm_particle_end", // output file, resulted particles after trmming
     "trm_boundary_end", // output file, resulted boundaries after trmming
     "dep_debug");       // output file, debug 
  */
  /*
  // degravitation, no boundary, quasi-static
  assemb.deGravitation(1000,               // total_steps
		       5,                  // number of snapshots
		       1,                  // print interval
		       true,               // recreate from input file
		       "dep_particle_end", // input file, initial particles
		       "dgr_particle",     // output file, resulted particles, including snapshots 
		       "dgr_contact",      // output file, resulted contacts, including snapshots 
		       "dgr_progress",     // output file, statistical info
		       "dgr_debug");       // output file, debug info
  // container properties
  dem::Rectangle container(0.05, 0.05, 0.05, vec(0, 0, 0)); // dimx, dimy, dimz, center
  // particle shape, size and percentage
  REAL ptcl_ratio_ba = 1;//0.8;  // ratio of radius b to radius a
  REAL ptcl_ratio_ca = 1;//0.6;  // ratio of radius c to radius a
  std::vector<REAL> percent; // mass percentage of particles smaller than a certain size
  std::vector<REAL> ptclSize;    // particle size
  percent.push_back(1.00); ptclSize.push_back(2.5e-3);
  //percent.push_back(0.80); ptclSize.push_back(2.3e-3);
  //percent.push_back(0.60); ptclSize.push_back(2.0e-3);
  //percent.push_back(0.30); ptclSize.push_back(1.5e-3);
  //percent.push_back(0.10); ptclSize.push_back(1.0e-3);
  dem::Gradation ptclGradation(percent.size(), percent, ptclSize, ptcl_ratio_ba, ptcl_ratio_ca);
  assemb.setContainer(container);
  assemb.setGradation(ptclGradation);
  ///////////////////////////////////////////////////////////////////////////////////
  // trim(), createMemParticle() and iso_MemBdry() must be called together because:
  // 1. trim() and createMemParticle() share variable HistoryNum.
  // 2. springs in iso_MemBdry() reference particles created in createMemParticle().
  assemb.trim(false,               // recreate from input file or not
	      "dgr_particle_end",  // input
	      "trm_particle_end"); // output
  assemb.createMemParticle(0.25,   // size relative to minimum radius
			   false,  // recreate from input file or not
			   "trm_particle_end",
			   "mem_particle_end");
  assemb.iso_MemBdry(100000,
		     100,
		     10,
		     1.0e+3,
		     0.25,        // size relative to minimum radius
		     false,       // recreate from input file or not, must be false
		     "mem_particle_end",
		     "iso_particle",
		     "iso_contact",
		     "iso_progress",
		     "iso_debug");
  */
  
  /*
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Part 3: record run time
  time(&time2);
  dem::g_timeinf << std::endl 
	    << "simulation start time: " << ctime(&time1);
  dem::g_timeinf << "simulation  end  time: " << ctime(&time2);
  dem::g_timeinf.close();
  */

  return 0;
}

 
////////////////////////////////////////////////////////////////////////////////////////////////////
// Notes:
//
// particleLayers (settings for free particles):
//  0 - a single free particle
//  1 - a layer of free particles
//  2 - multiple layers of free particles
//
// particle type (settings for individual particle):
//  0 - free particle
//  1 - fixed particle
//  2 - special case 2 (pure moment): translate first, then rotate only, MNT_START needs to be defined
//  3 - special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
//  4 - special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only
//  5 - free boundary particle
// 10 - ghost particle

////////////////////////////////////////////////////////////////////////////////////////////////////
// Various types of simulation, copy into part 2 of main() function to run, only ONE block a time.
/*

    // for triaxalPtclBdry
    assemb.TrimPtclBdryByHeight(0.061,
			   "dep_particle_end",
			   "dep_particle_trimmed");


    // for de-fe coupling			   
    assemb.TrimPtclBdryByHeight(0.061,
			   "pile_particle_ini",
			   "pile_particle_trimmed");


    assemb.triaxialPtclBdryIni(100000,             // total_steps
			  100,                // number of snapshots
                          10,                 // print interval
			  2.5e+6,             // confining pressure to achieve
			  "ini_particle_500k",// input file, initial particles
			  "ini_boundary_500k",// input file, initial boundaries
			  "ini_particle",     // output file, resulted particles, including snapshots 
			  "ini_boundary",     // output file, resulted boundaries
			  "ini_contact",      // output file, resulted contacts, including snapshots 
			  "ini_progress",     // output file, statistical info
			  "ini_debug");       // output file, debug info


    assemb.triaxialPtclBdryIni(100000,             // total_steps
			  100,                // number of snapshots
                          10,                 // print interval
			  5.0e+5,             // confining pressure to achieve
			  "ini_particle_ini", // input file, initial particles
			  "ini_boundary_ini", // input file, initial boundaries
			  "ini_particle",     // output file, resulted particles, including snapshots 
			  "ini_boundary",     // output file, resulted boundaries
			  "ini_contact",      // output file, resulted contacts, including snapshots 
			  "ini_progress",     // output file, statistical info
			  "ini_debug");       // output file, debug info


    assemb.triaxialPtclBdry(100000,             // total_steps
		       100,                // number of snapshots
                       10,                 // print interval
		       "tri_particle_ini", // input file, initial particles
		       "tri_boundary_ini", // input file, initial boundaries
		       "tri_particle",     // output file, resulted particles, including snapshots 
		       "tri_boundary",     // output file, resulted boundaries
		       "tri_contact",      // output file, resulted contacts, including snapshots 
		       "tri_progress",     // output file, statistical info
		       "tri_balanced",     // output file, balanced status
		       "tri_debug");       // output file, debug info


  // container properties
  dem::Rectangle container(0.05, 0.05, 0.05, vec(0, 0, 0)); // dimx, dimy, dimz, center
  // particle shape, size and percentage
  REAL ptcl_ratio_ba = 0.8;  // ratio of radius b to radius a
  REAL ptcl_ratio_ca = 0.6;  // ratio of radius c to radius a
  std::vector<REAL> percent; // mass percentage of particles smaller than a certain size
  std::vector<REAL> ptclSize;    // particle size
  percent.push_back(1.00); ptclSize.push_back(2.5e-3);
  //percent.push_back(0.80); ptclSize.push_back(2.3e-3);
  //percent.push_back(0.60); ptclSize.push_back(2.0e-3);
  //percent.push_back(0.30); ptclSize.push_back(1.5e-3);
  //percent.push_back(0.10); ptclSize.push_back(1.0e-3);
  dem::Gradation ptclGradation(percent.size(), percent, ptclSize, ptcl_ratio_ba, ptcl_ratio_ca);
  assemb.setContainer(container);
  assemb.setGradation(ptclGradation);
  assemb.depositIntoContainer(2,                  // particleLayers, setting of free particles 
			 100000,             // total_steps
			 100,                // number of snapshots
			 10,                 // print interval
			 3.0,                // relative height of floating particles based on container height
			 "flo_particle_end", // output file, initial particles for depositing
			 "dep_boundary_ini", // output file, initial boundaries for depositing
			 "dep_particle",     // output file, resulted particles, including snapshots 
			 "dep_contact",      // output file, resulted contacts, including snapshots 
			 "dep_progress",     // output file, statistical info
			 "trm_particle_end", // output file, resulted particles after trmming
			 "trm_boundary_end", // output file, resulted boundaries after trmming
			 "dep_debug");       // output file, debug 


    // size, shape, and Gradation of particles
    int rorc             = 1;     // rectangular = 1 or cylindrical = 0
    REAL dimn     = 0.05;  // specimen dimension
    REAL ratio_ba = 0.8;   // ratio of radius b to radius a
    REAL ratio_ca = 0.6;   // ratio of radius c to radius a
    std::vector<REAL> percent;  // mass percentage of particles smaller than a certain size
    std::vector<REAL> ptclSize; // particle size
    percent.push_back(1.00); ptclSize.push_back(2.5e-3);
    //percent.push_back(0.80); ptclSize.push_back(2.0e-3);
    //percent.push_back(0.60); ptclSize.push_back(1.6e-3);
    //percent.push_back(0.30); ptclSize.push_back(1.0e-3);
    //percent.push_back(0.10); ptclSize.push_back(0.5e-3);
    dem::Gradation grad(rorc, dimn, ratio_ba, ratio_ca, percent.size(), percent, ptclSize);
    assemb.deposit_PtclBdry(grad,
		       2,                  // particleLayers, setting of free particles
		       1.0,                // relative container size, 0.8/1.0/1.2---small/medium/large
		       100000,             // total_steps
		       100,                // number of snapshots
                       10,                 // print interval
		       "flo_particle_end", // output file, initial particles
		       "dep_particle",     // output file, resulted particles, including snapshots 
		       "dep_contact",      // output file, resulted contacts, including snapshots 
		       "dep_progress",     // output file, statistical info
		       "dep_debug");       // output file, debug info
 
   
    assemb.scale_PtclBdry(20000,             // total_steps
		     100,               // number of snapshots  
                     10,                // print interval
		     0.05,              // dimension of particle-composed-boundary
		     1.0,               // relative container size, 0.8/1.0/1.2---small/medium/large
		     "dep_particle_end",// input file, initial particles
		     "scl_particle",    // output file, resulted particles, including snapshots 
		     "scl_contact",     // output file, resulted contacts, including snapshots
		     "scl_progress",    // output file, statistical info
		     "scl_debug");      // output file, debug info


    assemb.ellipPile_Disp(50000,              // total_steps
		     100,                // number of snapshots
                     10,                 // print interval
		     0.05,               // dimension of particle-composed-boundary
		     1.0,                // relative container size, 0.8/1.0/1.2---small/medium/large
		     "pile_particle_ini",// input file, initial particles, an ellipsoidal pile info added
		     "pile_particle",    // output file, resulted particles, including snapshots 
		     "pile_contact",     // output file, resulted contacts, including snapshots 
		     "pile_progress",    // output file, statistical info
		     "pile_debug");      // output file, debug info


    assemb.rectPile_Disp(50000,              // total_steps
		    100,                // number of snapshots
                    10,                 // print interval
		    "pile_particle_ini",// input file, initial particles
		    "pile_boundary_ini",// input file, initial boundaries, rectangular pile boundary info added
		    "pile_particle",    // output file, resulted particles, including snapshots 
		    "pile_boundary",    // output file, resulted boundaries
		    "pile_contact",     // output file, resulted contacts, including snapshots 
		    "pile_progress",    // output file, statistical info
		    "pile_debug");      // output file, debug info


    assemb.ellipPile_Impact(50000,              // total_steps
		       100,                // number of snapshots
                       10,                 // print interval
		       0.05,               // size of particle-composed-boundary
		       "ipt_particle_ini", // input file, initial particles, an ellipsoidal pile info added
		       "dep_boundary_ini", // input file, initial boundaries
		       "ipt_particle",     // output file, resulted particles, including snapshots 
		       "ipt_contact",      // output file, resulted contacts, including snapshots 
		       "ipt_progress",     // output file, statistical info
		       "ipt_debug");       // output file, debug info


    assemb.ellipPile_Impact_p(50000,              // total_steps
			 100,                // number of snapshots
                         10,                 // print interval
			 0.05,               // size of particle-composed-boundary
			 "ipt_particle_ini", // input file, initial particles, an ellipsoidal pile info added
			 "ipt_particle",     // output file, resulted particles, including snapshots 
			 "ipt_contact",      // output file, resulted contacts, including snapshots 
			 "ipt_progress",     // output file, statistical info
			 "ipt_debug");       // output file, debug info


  assemb.deposit(100000,             // total_steps
		 100,                // number of snapshots
		 10,                 // print interval
		 "flo_particle_end", // input file, initial particles
		 "dep_boundary_ini", // input file, initial boundaries
		 "dep_particle",     // output file, resulted particles, including snapshots 
		 "dep_contact",      // output file, resulted contacts, including snapshots 
		 "dep_progress",     // output file, statistical info
		 "dep_debug");       // output file, debug info    


    assemb.squeeze(300000,             // total_steps
	      100000,             // initial_steps to reach equilibrium
              100,                // number of snapshots
              10,                 // print interval
	      +1,                 // -1 squeeze; +1 loosen
	      "flo_particle_end", // input file, initial particles
	      "dep_boundary_ini", // input file, initial boundaries
	      "dep_particle",     // output file, resulted particles, including snapshots 
	      "dep_boundary",     // output file, resulted boundaries
	      "dep_contact",      // output file, resulted contacts, including snapshots 
	      "dep_progress",     // output file, statistical info
	      "dep_debug");       // output file, debug info


    assemb.collapse(rorc,
	       100000,
	       100,
               10,                // print interval
	       "cre_particle",    // input file, initial particles
	       "clp_boundary",    // output file, initial boundaries
	       "clp_particle",    // output file, resulted particles, including snapshots
	       "clp_contact",     // output file, resulted contacts, including snapshots 
	       "clp_progress",    // output file, statistical info
	       "clp_debug");      // output file, debug info


    assemb.isotropic(100000,             // total_steps
		100,                // number of snapshots
                10,                 // print interval
		1.0e+3,             // a low confining pressure to achieve for initialization
		"cre_particle",     // input file, initial particles
		"cre_boundary",     // input file, initial boundaries
		"iso_particle",     // output file, resulted particles, including snapshots 
		"iso_boundary",     // output file, resulted boundaries 
		"iso_contact",      // output file, resulted contacts, including snapshots 
		"iso_progress",     // output file, statistical info
		"iso_balanced",     // output file, balanced status
		"iso_debug");       // output file, debug info

    
    assemb.isotropic(100000,             // total_steps
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

    
    assemb.isotropic(100000,             // total_steps
		100,                // number of snapshots
                10,                 // print interval
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

    
    REAL sigma_values[4]={1.0e+5, 5.0e+5, 1.0e+5, 7.0e+5}; // last one must be a larger value	
    assemb.isotropic(100000,             // total_steps
		100,                // number of snapshots
                10,                 // print interval
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

    
    assemb.odometer(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // print interval
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

    
    REAL sigma_values[4]={1.0e+5, 5.0e+5, 1.0e+5, 1.0e+6}; // last one must be a larger value	
    assemb.odometer(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // print interval
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


    assemb.triaxial(100000,             // total_steps
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

    
    assemb.triaxial(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // print interval
	       3.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_300k",// input file, initial particles
	       "iso_boundary_300k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary"      // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    assemb.triaxial(100000,             // total_steps
	       100,                // number of snapshots
               10,                 // print interval
	       5.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_500k",// input file, initial particles
	       "iso_boundary_500k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary",     // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    assemb.triaxial(120000,             // total_steps
	       30000,              // at which step to unload
	       100,                // number of snapshots
               10,                 // print interval
	       3.0e+5,             // pre-existing confining pressure from initial isotropic compression
	       "iso_particle_300k",// input file, initial particles
	       "iso_boundary_300k",// input file, initial boundaries
	       "tri_particle",     // output file, resulted particles, including snapshots 
	       "tri_boundary",     // output file, resulted boundaries
	       "tri_contact",      // output file, resulted contacts, including snapshots 
	       "tri_progress",     // output file, statistical info
	       "tri_balanced",     // output file, balanced status
	       "tri_debug");       // output file, debug info

    
    assemb.truetriaxial(100000,             // total_steps
		   100,                // number of snapshots
                   10,                 // print interval
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


    assemb.unconfined(100000,             // total_steps
		 100,                // number of snapshots
                 10,                 // print interval
		 "flo_particle_end", // input file, initial particles
		 "unc_boundary",     // input file, initial boundaries
		 "unc_particle",     // output file, resulted particles, including snapshots 
		 "unc_contact",      // output file, resulted contacts, including snapshots 
		 "unc_progress",     // output file, statistical info
		 "unc_debug");       // output file, debug info
*/
