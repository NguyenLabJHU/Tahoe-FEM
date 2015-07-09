#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "realtypes.h"
#include "Parameter.h"
#include "Vec.h"
#include "Gradation.h"
#include "Particle.h"
#include "Contact.h"
#include "Boundary.h"
#include "Particle.h"
#include "Rectangle.h"
#include "Cylinder.h"
#include "Spring.h"
#include "Fluid.h"
#include <cstddef>
#include <map>
#include <vector>
#include <fstream>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace dem {
  
  class Assembly {  
  private:

    // particles property
    Gradation               gradation;       // particles gradation, broadcast among processes for once
    std::vector<Particle *> allParticleVec;  // all particles, only meaningful to root process
    std::vector<Particle *> particleVec;     // particles per process
    std::size_t             trimHistoryNum;  // historical maximum numbering before trimming, only meaningful to root process

    std::vector<Contact>    contactVec;      // contacts per process
    std::vector<ContactTgt> contactTgtVec;   // tangential contact force and displacement per process
    std::size_t             allContactNum;   // estimated total contact number, only meaningful to root process
    
    std::vector< std::vector< std::vector<Particle *> > > memBoundary; // membrane particle boundaries
    std::vector<Spring *>   springVec;       // springs connecting membrane particles

    // container property
    Rectangle allContainer;// whole container, broadcast among processes for once
    Rectangle container;   // container per process
    Rectangle cavity;      // cavity inside container
    Rectangle grid;        // adaptive compute grid, broadcast among processes for once, updated per process
    
    // boundaries property
    std::vector<Boundary *> boundaryVec;       // rigid boundaries, broadcast among processes upon changed.
    std::vector<Boundary *> mergeBoundaryVec;  // rigid boundaries with stats from all processes
    std::vector<Boundary *> cavityBoundaryVec; // rigid cavity boundaries
    std::map<std::size_t,std::vector<BoundaryTgt> > boundaryTgtMap; // particle-boundary contact tangential info
   
    // fluid property
    Fluid fluid;

    // average data
    REAL avgNormal;        // only meaningful to root process
    REAL avgShear;         // only meaningful to root process
    REAL avgPenetr;        // only meaningful to root process

    // energy data
    REAL transEnergy;      // only meaningful to root process
    REAL rotatEnergy;      // only meaningful to root process
    REAL kinetEnergy;      // only meaningful to root process
    REAL graviEnergy;      // only meaningful to root process
    REAL mechaEnergy;      // only meaningful to root process

    // time step
    REAL vibraTimeStep;    // meaningful to all processes
    REAL impactTimeStep;   // meaningful to all processes

    // MPI data
    boost::mpi::communicator boostWorld;
    MPI_Comm mpiWorld, cartComm;
    std::vector<std::size_t> bdryProcess;
    int mpiProcX, mpiProcY, mpiProcZ;
    int mpiRank, mpiSize, mpiTag, mpiCoords[3];
    int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
    int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2; 
    int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2; 
    int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2; 
    int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2; 
    int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;
    std::vector<Particle *> rParticleX1, rParticleX2; // r stands for received
    std::vector<Particle *> rParticleY1, rParticleY2; 
    std::vector<Particle *> rParticleZ1, rParticleZ2; 
    std::vector<Particle *> rParticleX1Y1, rParticleX1Y2, rParticleX1Z1, rParticleX1Z2; 
    std::vector<Particle *> rParticleX2Y1, rParticleX2Y2, rParticleX2Z1, rParticleX2Z2; 
    std::vector<Particle *> rParticleY1Z1, rParticleY1Z2, rParticleY2Z1, rParticleY2Z2; 
    std::vector<Particle *> rParticleX1Y1Z1, rParticleX1Y1Z2, rParticleX1Y2Z1, rParticleX1Y2Z2; 
    std::vector<Particle *> rParticleX2Y1Z1, rParticleX2Y1Z2, rParticleX2Y2Z1, rParticleX2Y2Z2; 
    std::vector<Particle *> recvParticleVec;  // received particles per process
    std::vector<Particle *> mergeParticleVec; // merged particles per process
      
    // stream
    std::ofstream progressInf;
    std::ofstream balancedInf;

  public:
    Assembly()
      :trimHistoryNum(0), allContactNum(0), avgNormal(0), avgShear(0), avgPenetr(0),
      transEnergy(0), rotatEnergy(0), kinetEnergy(0), graviEnergy(0), mechaEnergy(0), 
      vibraTimeStep(0), impactTimeStep(0)
      {}
    
    ~Assembly() {
      // release memory pointed to by pointers in the container
      for(std::vector<Particle *>::iterator it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
	delete (*it);

      for(std::vector<Particle *>::iterator it = particleVec.begin(); it != particleVec.end(); ++it)
	delete (*it);

      for(std::vector<Boundary *>::iterator it = boundaryVec.begin(); it != boundaryVec.end(); ++it)
	delete (*it);

      for(std::vector<Boundary *>::iterator it = cavityBoundaryVec.begin(); it != cavityBoundaryVec.end(); ++it)
	delete (*it);

      for(std::vector<Spring *>::iterator it = springVec.begin(); it != springVec.end(); ++it)
	delete (*it); 

      // in case of consecutive simulations
      allParticleVec.clear();
      particleVec.clear();
      boundaryVec.clear();
      cavityBoundaryVec.clear();
      springVec.clear();
    }
   
    void setCommunicator(boost::mpi::communicator &comm);
    void setContainer(Rectangle cont) { allContainer = cont; } 
    void setGrid(Rectangle cont) { grid = cont; } 
    void setGradation(Gradation grad) { gradation = grad; }

    void tuneMassPercent();	// July 3, 2015
    void calcMassPercent();	// July 3, 2015
    void depositIntoContainer(); 
    void resumeDepositIntoContainer();
    void expandCavityParticle();
    void resumeExpandCavityParticle();
    void generateParticle(std::size_t particleLayers,
			  const char *genParticle);
    void buildBoundary(std::size_t boundaryNum,
		       const char *boundaryFile);
    void trimOnly();
    void trim(bool toRebuild,
	      const char *inputParticle,
	      const char *trmParticle);
    void removeBySphere();
    void deposit(const char *inputBoundary,
		 const char *inputParticle);
    void proceedFromPreset();
    void coupleWithGas();    

    void isotropic();
    void odometer();
    void triaxial();
    void planeStrain();
    void trueTriaxial();
    bool tractionErrorTol(REAL sigma, std::string type, REAL sigmaX=0, REAL sigmaY=0);
    void getStartDimension(REAL &distX, REAL &distY, REAL &distZ);

    void setCavity(Rectangle cav) { cavity = cav; }

    void readParticle(const char *str);	// July 3, 2015
    void readBoundary(const char *str);
    void scatterParticle();
    void commuParticle();
    bool isBdryProcess();
    void releaseRecvParticle();
    void releaseGatheredParticle();
    void releaseGatheredContact();
    void migrateParticle();
    void removeParticleOutRectangle();
    void gatherParticle();
    void gatherBdryContact();

    void updateGrid();
    void updateGridMinX();
    void updateGridMaxX();
    void updateGridMinY();
    void updateGridMaxY();
    void updateGridMinZ();
    void updateGridMaxZ();    

    void openDepositProg(std::ofstream &ofs, const char *str);
    void printDepositProg(std::ofstream &ofs);
    void openCompressProg(std::ofstream &ofs, const char *str);
    void printCompressProg(std::ofstream &ofs, REAL distX, REAL distY, REAL distZ);
    void openParticleProg(std::ofstream &ofs, const char *str);
    void closeProg(std::ofstream &ofs);

    void trimCavity(bool toRebuild, const char *Particlefile, const char *cavParticle);
    void readCavityBoundary(const char *boundaryfile);
    void buildCavityBoundary(std::size_t existMaxId, const char *boundaryfile);
    void findContact();                           // detect and resolve contact between particles, July 3, 2015
    void findBdryContact();                       // find particles on boundaries
    void findParticleOnCavity();                  // find particle on cavity boundaries
    
    void clearContactForce();                     // clear forces and moments for all particles
    void internalForce();                         // calculate inter-particle forces
    void dragForce();
    void springForce();
    void boundaryForce();                         // calcualte forces between rigid boundaries and particles
    void cavityBoundaryForce();
    void updateParticle();                        // update motion of particles
    
    REAL ellipPileForce();                        // for force pile only
    void ellipPileUpdate();                       // for force pile only
    
    Vec  ellipPileDimn();
    REAL ellipPileTipZ();
    REAL ellipPilePeneVol();
  
    void updateBoundary(REAL simga, std::string type, REAL sigmaX=0, REAL sigmaY=0);
    
    REAL getMass() const; 
    REAL getAvgPenetr() const;
    REAL getParticleVolume() const;

    void calcTimeStep();
    void calcVibraTimeStep();
    void calcImpactTimeStep();
    void calcContactNum();

    REAL getAvgTransVelocity() const;
    REAL getAvgRotatVelocity() const;
    REAL getAvgForce() const;
    REAL getAvgMoment() const;

    void calcTransEnergy();
    void calcRotatEnergy();
    void calcKinetEnergy();
    void calcGraviEnergy(REAL ref);
    void calcMechaEnergy();
    void gatherEnergy();
    
    void setTrimHistoryNum(std::size_t n) { trimHistoryNum = n; }
    void printParticle(const char *str) const; // print all particles, July 3, 2015
    void printBdryContact(const char *str) const; // print all boundary contact info
    void printParticle(const char *str, std::vector<Particle *>  &particleVec) const; // print particles info, July 3, 2015
    void printMemParticle(const char *str) const; // print membrane particles
    void plotSpring(const char *str) const;    // print springs in Tecplot format
    void plotBoundary(const char *str) const;
    void plotGrid(const char *str) const;
    void plotCavity(const char *str) const;
    void checkMembrane(std::vector<REAL> &vx ) const;
    void printContact(char *str) const;        // print contacts information
    void printBoundary(const char *str) const; // print rigid boundaries info
    void printCavityBoundary(const char *str) const; // print cavity boundaries
    void printCavityParticle(std::size_t total, const char *str) const;
    
  // continue to deposit after a cavity is created inside the particle assemblage
  void depositAfterCavity(std::size_t  total_steps,  
			  std::size_t  snapNum,
			  std::size_t  interval,
			  const char *iniptclfile,   
			  const char *inibdryfile,
			  const char *inicavefile,
			  const char *Particlefile, 
			  const char *contactfile,
			  const char *progressfile, 
			  const char *debugfile);

  // create a specimen by depositing particles into particle boundaries
  void deposit_PtclBdry(Gradation &grad,
			std::size_t  freetype,
			REAL  rsize,
			std::size_t  total_steps,  
			std::size_t  snapNum,
			std::size_t  interval,
			const char *iniptclfile,   
			const char *Particlefile, 
			const char *contactfile,
			const char *progressfile, 
			const char *debugfile);
  
  // scale the assembly with particle boundaries from deposited state until it reaches steady state
  void scale_PtclBdry(std::size_t        total_steps  =50000,             // total_steps
		      std::size_t        snapNum    =100,               // number of snapNum   
		      std::size_t        interval     =10,                // print interval
		      REAL        dimn         =0.05,              // dimension of particle-composed-boundary
		      REAL        rsize        =1.0,               // relative container size
		      const char *iniptclfile  ="dep_particle_end",// input file, initial particles
		      const char *Particlefile ="scl_particle",    // output file, resulted particles, including snapNum 
		      const char *contactfile  ="scl_contact",     // output file, resulted contacts, including snapNum
		      const char *progressfile ="scl_progress",    // output file, statistical info
		      const char *debugfile    ="scl_debug");      // output file, debug info
  
  
  // generate particles in space for particle boundaries
  void generate_p(Gradation &grad,
		  const char *str,
		  std::size_t freetype,
		  REAL rsize,
		  REAL ht);
  
 
  void deGravitation(std::size_t  total_steps,  
		     std::size_t  snapNum,
		     std::size_t  interval,
		     bool  toRebuild,
		     const char *iniptclfile,   
		     const char *Particlefile, 
		     const char *contactfile,
		     const char *progressfile, 
		     const char *debugfile);
  
  // actual deposit function for particle boundaries
  void deposit_p(std::size_t        total_steps  =50000,             // total_steps
		 std::size_t        snapNum    =100,               // number of snapNum   
		 std::size_t        interval     =10,                // print interval 
		 REAL dimn   =0.05,                           // dimension of particle-composed-boundary
		 REAL rsize  =1.0,                            // relative container size
		 const char *iniptclfile  ="flo_particle_end",// input file, initial particles
		 const char *Particlefile ="dep_particle",    // output file, resulted particles, including snapNum 
		 const char *contactfile  ="dep_contact",     // output file, resulted contacts, including snapNum
		 const char *progressfile ="dep_progress",    // output file, statistical info
		 const char *debugfile    ="dep_debug");      // output file, debug info
  
  //squeeze paticles inside a container by moving the boundaries
  void squeeze(std::size_t        total_steps  =20000,               // total_steps
	       std::size_t        init_steps   =5000,                // initial_steps to reach equilibrium
	       std::size_t        snapNum      =100,                 // number of snapNum   
	       std::size_t        interval     =10,                  // print interval 
	       int                flag         =-1,                  // -1 squeeze; +1 loosen
	       const char *iniptclfile  ="flo_particle_end",  // input file, initial particles
	       const char *inibdryfile  ="dep_boundary_ini",  // input file, initial boundaries
	       const char *Particlefile ="dep_particle",      // output file, resulted particles, including snapNum 
	       const char *boundaryfile ="dep_boundary",      // output file, resulted boundaries
	       const char *contactfile  ="dep_contact",       // output file, resulted contacts, including snapNum
	       const char *progressfile ="dep_progress",      // output file, statistical info
	       const char *debugfile    ="dep_debug");        // output file, debug info
  
  void deposit_repose(std::size_t  interval,
		      const char *inibdryfile,
		      const char *Particlefile, 
		      const char *contactfile,
		      const char *progressfile, 
		      const char *debugfile);
  
  void angleOfRepose(std::size_t  interval,
		     const char *inibdryfile,
		     const char *Particlefile, 
		     const char *contactfile,
		     const char *progressfile, 
		     const char *debugfile);
  
  REAL getPtclMinX(const std::vector<Particle *> &particleVec) const;
  REAL getPtclMaxX(const std::vector<Particle *> &particleVec) const;
  REAL getPtclMinY(const std::vector<Particle *> &particleVec) const;
  REAL getPtclMaxY(const std::vector<Particle *> &particleVec) const;
  REAL getPtclMinZ(const std::vector<Particle *> &particleVec) const;
  REAL getPtclMaxZ(const std::vector<Particle *> &particleVec) const;
  
  void collapse(std::size_t  total_steps,  
		std::size_t  snapNum,
		std::size_t  interval,
		const char *iniptclfile,
		const char *initboundary,
		const char *Particlefile,
		const char *contactfile,
		const char *progressfile,
		const char *debugfile);
  
  void createMemParticle(REAL rRadius,
			 bool toRebuild,
			 const char *Particlefile,
			 const char *allParticle);
  
  void iso_MemBdry(std::size_t  total_steps,  
		   std::size_t  snapNum, 
		   std::size_t  interval,
		   REAL  sigma3,
		   REAL  rRadius,
		   bool  toRebuild,
		   const char *iniptclfile, 
		   const char *Particlefile,
		   const char *contactfile, 
		   const char *progressfile,
		   const char *debugfile);
  
  void TrimPtclBdryByHeight(REAL height,
			    const char *iniptclfile,
			    const char *Particlefile);
  
  void applyParticleBoundary(std::size_t         total_steps  =100000,
			     std::size_t         snapNum    =100,
			     std::size_t         nterval      =10,
			     REAL         sigma        =1.0e+4,
			     const char *iniptclfile  ="cre_particle",
			     const char *inibdryfile  ="cre_bounary",
			     const char *Particlefile ="iso_particle",
			     const char *boundaryfile ="iso_boundary",
			     const char *contactfile  ="iso_contact",
			     const char *progressfile ="iso_progress",
			     const char *balancedfile ="iso_balanced",
			     const char *debugfile    ="iso_debug");
  
  // The confining pressure is 500kPa. This function initializes triaxial compression test.
  void triaxialPtclBdryIni(std::size_t         total_steps  =10000,
			   std::size_t         snapNum    =100,
			   std::size_t         interval     =10,
			   REAL         sigma        =5.0e+5,
			   const char *iniptclfile  ="ini_particle_ini",
			   const char *inibdryfile  ="ini_boundary_ini",
			   const char *Particlefile ="ini_particle", 
			   const char *boundaryfile ="ini_boundary", 
			   const char *contactfile  ="ini_contact",
			   const char *progressfile ="ini_progress",
			   const char *debugfile    ="ini_debug");
  
  // The confining pressure is 500kPa. This function performs triaxial compression test.
  // Displacement boundaries are used in axial direction.
  void triaxialPtclBdry(std::size_t         total_steps  =100000,
			std::size_t         snapNum    =100,
			std::size_t         interval     =10,
			const char *iniptclfile  ="iso_particle_100k",
			const char *inibdryfile  ="iso_boundary_100k",
			const char *Particlefile ="tri_particle", 
			const char *boundaryfile ="tri_boundary", 
			const char *contactfile  ="tri_contact",
			const char *progressfile ="tri_progress",
			const char *balancedfile ="tri_balanced", 
			const char *debugfile    ="tri_debug");
   
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // A rectangular pile is then drived into the particles using displacement control.
  void rectPile_Disp(std::size_t         total_steps  =50000,
		     std::size_t         snapNum    =100,
		     std::size_t         interval     =10,
		     const char *iniptclfile  ="pile_particle_ini",
		     const char *inibdryfile  ="pile_boundary_ini",
		     const char *Particlefile ="pile_particle", 
		     const char *boundaryfile ="pile_boundary", 
		     const char *contactfile  ="pile_contact",
		     const char *progressfile ="pile_progress",
		     const char *debugfile    ="pile_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // An ellipsoidal pile is then drived into the particles using displacement control.
  void ellipPile_Disp(std::size_t        total_steps  =50000,  
		      std::size_t        snapNum    =100, 
		      std::size_t         interval     =10,
		      REAL dimn         =0.05,
		      REAL rsize        =1.0,
		      const char *iniptclfile  ="pile_particle_ini",
		      const char *Particlefile ="pile_particle", 
		      const char *contactfile  ="pile_contact",  
		      const char *progressfile ="pile_progress",
		      const char *debugfile    ="pile_debug");
  
  // The specimen has been deposited with gravitation within rigid boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
  void ellipPile_Impact(std::size_t        total_steps  =50000,  
			std::size_t        snapNum    =100, 
			std::size_t        interval     =10,
			REAL dimn         =0.05,
			const char *iniptclfile  ="ipt_particle_ini",
			const char *inibdryfile  ="dep_boundary_ini",
			const char *Particlefile ="ipt_particle", 
			const char *contactfile  ="ipt_contact",  
			const char *progressfile ="ipt_progress",
			const char *debugfile    ="ipt_debug");
  
  // The specimen has been deposited with gravitation within particle boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
  void ellipPile_Impact_p(std::size_t        total_steps  =50000,  
			  std::size_t        snapNum    =100, 
			  std::size_t        interval     =10,
			  REAL dimn         =0.05,
			  const char *iniptclfile  ="ipt_particle_ini",
			  const char *Particlefile ="ipt_particle", 
			  const char *contactfile  ="ipt_contact",  
			  const char *progressfile ="ipt_progress",
			  const char *debugfile    ="ipt_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // An ellipsoidal pile is then drived into the particles using force control.
  void ellipPile_Force(std::size_t        total_steps  =50000,  
		       std::size_t        snapNum    =100, 
		       std::size_t        interval     =10,
		       REAL dimn         =0.05,
		       REAL force        =1.0e+4,
		       std::size_t  division           =100,
		       const char *iniptclfile  ="pile_particle_ini",
		       const char *Particlefile ="pile_particle", 
		       const char *contactfile  ="pile_contact",  
		       const char *progressfile ="pile_progress",
		       const char *balancedfile ="pile_balanced",
		       const char *debugfile    ="pile_debug");

  public:
    void findParticleInRectangle(const Rectangle &container,
				 const std::vector<Particle *> &allParticle,
				 std::vector<Particle *> &foundParticle);
    
  };
  
} // namespace dem

#endif
