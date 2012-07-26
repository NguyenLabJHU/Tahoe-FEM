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
#include <map>
#include <vector>
#include <fstream>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace dem {
  
  class Assembly {
    typedef Contact<Particle>  CONTACT;
    typedef Boundary<Particle> BOUNDARY;
    
  private:

    // particles property
    int  trimHistoryNum;                     // historical maximum numbering before trimming
    int  possContactNum;                     // possible contact number based on spherical distances
    int  actualContactNum;                   // actual contact number based on contact resolution
    Gradation               gradation;       // particles gradation
    std::vector<Particle*>  allParticleVec;  // all particles
    std::vector<Particle*>  particleVec;     // particles per process
    
    std::vector<CONTACT>    contactVec;      // contacts per process
    std::vector<ContactTgt> allContactTgtVec;// all tangential contact force and displacement
    std::vector<ContactTgt> contactTgtVec;   // tangential contact force and displacement per process
    
    std::vector< std::vector< std::vector<Particle*> > > memBoundary; // membrane particle boundaries
    std::vector<Spring*>    springVec;       // springs connecting membrane particles
    
    // container property
    Rectangle allContainer;// whole container
    Rectangle container;   // container per process
    Rectangle cavity;      // cavity inside container
    
    // boundaries property
    std::vector<BOUNDARY*> boundaryVec;       // rigid boundaries
    std::vector<BOUNDARY*> cavityBoundaryVec; // rigid cavity boundaries
    std::map<int,std::vector<BoundaryTgt> > boundaryTgtMap; // particle-boundary contact tangential info
    
    // MPI data
    boost::mpi::communicator boostWorld;
    MPI_Comm mpiWorld, cartComm;
    int mpiProcX, mpiProcY, mpiProcZ;
    int mpiRank, mpiSize, mpiTag, mpiCoords[3];
    int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
    int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2; 
    int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2; 
    int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2; 
    int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2; 
    int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;
    std::vector<Particle*> rParticleX1, rParticleX2; // r stands for received
    std::vector<Particle*> rParticleY1, rParticleY2; 
    std::vector<Particle*> rParticleZ1, rParticleZ2; 
    std::vector<Particle*> rParticleX1Y1, rParticleX1Y2, rParticleX1Z1, rParticleX1Z2; 
    std::vector<Particle*> rParticleX2Y1, rParticleX2Y2, rParticleX2Z1, rParticleX2Z2; 
    std::vector<Particle*> rParticleY1Z1, rParticleY1Z2, rParticleY2Z1, rParticleY2Z2; 
    std::vector<Particle*> rParticleX1Y1Z1, rParticleX1Y1Z2, rParticleX1Y2Z1, rParticleX1Y2Z2; 
    std::vector<Particle*> rParticleX2Y1Z1, rParticleX2Y1Z2, rParticleX2Y2Z1, rParticleX2Y2Z2; 
    std::vector<Particle*> recvParticleVec;  // received particles per process
    std::vector<Particle*> mergeParticleVec; // merged particles per process
    
    
  public:
    Assembly()
      :trimHistoryNum(0),
      possContactNum(0),
      actualContactNum(0)
      {}
    
    ~Assembly() {
      // release memory pointed to by pointers in the container
      std::vector<Particle*>::iterator pt;
      std::vector<BOUNDARY*>::iterator rt;
      std::vector<Spring*>::iterator   st;

      for(pt = allParticleVec.begin(); pt != allParticleVec.end(); ++pt)
	delete (*pt);

      for(pt = particleVec.begin(); pt != particleVec.end(); ++pt)
	delete (*pt);

      for(rt = boundaryVec.begin(); rt != boundaryVec.end(); ++rt)
	delete (*rt);

      for(rt = cavityBoundaryVec.begin(); rt != cavityBoundaryVec.end(); ++rt)
	delete (*rt);

      for(st = springVec.begin(); st != springVec.end(); ++st)
	delete (*st);    

      // in case of consecutive simulations
      allParticleVec.clear();
      particleVec.clear();
      boundaryVec.clear();
      cavityBoundaryVec.clear();
      springVec.clear();
    }
   
    void setCommunicator(boost::mpi::communicator &comm);
    void setContainer(Rectangle cont) { allContainer = cont; } 
    void setGradation(Gradation grad) { gradation = grad; }

    void depositIntoContainer(); 
    void generateParticle(int particleLayers,
			  const char *genParticle);
    void buildBoundary(int boundaryNum,
		       const char* boundaryFile);
    void buildBoundary(const char* boundaryFile);
    void trim(bool toRebuild,
	      const char* inputParticle,
	      const char* trmParticle);
    void deposit(int totalSteps,  
		 int snapNum,
		 int statInterv,
		 const char *inputBoundary,
		 const char *inputParticle);
    
    void setCavity(Rectangle cavi) {cavity = cavi;}

    void readParticle(const char* str);
    void readBoundary(const char* str);
    void scatterParticle();
    void commuParticle();
    bool isBdryProcess();
    void releaseRecvParticle();
    void transferParticle();
    void removeParticleOutRectangle();
    void gatherParticle();
    void updateContainerMinX();
    void updateContainerMaxX();
    void updateContainerMinY();
    void updateContainerMaxY();
    void updateContainerMinZ();
    void updateContainerMaxZ();    

    void trimCavity(bool toRebuild, const char* Particlefile, const char* cavParticle);
    void readCavityBoundary(const char* boundaryfile);
    void buildCavityBoundary(int existMaxId, const char* boundaryfile);
    void findContact();                           // detect and resolve contact between particles
    void findBdryContact();                       // find particles on boundaries
    void findParticleOnCavity();                  // find particle on cavity boundaries
    
    void clearContactForce();                     // clear forces and moments for all particles
    void internalForce();                         // calculate inter-particle forces
    void springForce();
    void boundaryForce();                         // calcualte forces between rigid boundaries and particles
    void boundaryForce(REAL penetr[],int cntnum[]);
    void cavityBoundaryForce();
    void updateParticle();                        // update motion of particles
    
    REAL ellipPileForce();                        // for force pile only
    void ellipPileUpdate();                       // for force pile only
    
    Vec  ellipPileDimn();
    REAL ellipPileTipZ();
    REAL ellipPilePeneVol();
  
    // if bn[i]=2, the 2nd rigid boundary should be updated according to rbctl[i],
    // totally num rigid boundaries must be updated
    void updateBoundary(int bn[], UPDATECTL rbctl[], int num);
    void updateBoundary6();
    
    REAL getDensity() const; 
    int  getPossContactNum() const {return  possContactNum;};
    int  getActualContactNum() const {return actualContactNum;}
    REAL getAvgPenetration() const;
    REAL getVibraTimeStep() const;
    REAL getImpactTimeStep() const;
    REAL getAvgVelocity() const;
    REAL getAvgForce() const;
    REAL getAvgOmga() const;
    REAL getAvgMoment() const;
    REAL getParticleVolume() const;
    Vec  getTopFreeParticlePosition() const;
    REAL getTransEnergy() const;
    REAL getRotatEnergy() const;
    REAL getKinetEnergy() const;
    REAL getPotenEnergy(REAL ref) const;
    
    Vec  getNormalForce(int bdry) const;       // get normal force acting on the bdry_th rigid boundary
    Vec  getShearForce(int bdry) const;        // get shear force acting on the bdry_th rigid boundary
    REAL getAvgNormal(int bdry) const;
    Vec  getApt(int bdry) const;               // get a point on bdry_th rigid boundary
    Vec  getDirc(int bdry) const;              // get the dirc of bdry_th rigid boundry
    REAL getArea(int bdry) const;
    REAL getAvgPressure() const;
    void setArea(int bdry,REAL a);             // set the area of the bdry-th rigid boundary be a
    void setTrimHistoryNum(int n) {trimHistoryNum = n;}
    void printParticle(const char* str) const; // print all particles info into a disk file
    void printParticle(const char* str, std::vector<Particle*>  &particleVec) const; // print particles info into a disk file
    void printMemParticle(const char* str) const; // print membrane particles info into a disk file
    void plotSpring(const char *str) const;    // print springs in Tecplot format
    void plotBoundary(const char *str) const;
    void plotGrid(const char *str) const;
    void plotCavity(const char *str) const;
    void checkMembrane(std::vector<REAL> &vx ) const;
    void printContact(const char* str) const;  // print contacts information
    void printBoundary(const char* str) const; // print rigid boundaries info to a disk file
    void printCavityBoundary(const char* str) const; // print cavity boundaries
    void printCavityParticle(int total, const char* str) const;
    
    void expandCavityParticles(bool toRebuild,
			       REAL percent,
			       const char* cavityptclfile,
			       const char* Particlefile,
			       const char* newptclfile);
    
  // continue to deposit after a cavity is created inside the particle assemblage
  void depositAfterCavity(int   total_steps,  
			  int   snapNum,
			  int   interval,
			  const char* iniptclfile,   
			  const char* inibdryfile,
			  const char* inicavefile,
			  const char* Particlefile, 
			  const char* contactfile,
			  const char* progressfile, 
			  const char* debugfile);

  // create a specimen by depositing particles into particle boundaries
  void deposit_PtclBdry(Gradation& grad,
			int   freetype,
			REAL  rsize,
			int   total_steps,  
			int   snapNum,
			int   interval,
			const char* iniptclfile,   
			const char* Particlefile, 
			const char* contactfile,
			const char* progressfile, 
			const char* debugfile);
  
  // scale the assembly with particle boundaries from deposited state until it reaches steady state
  void scale_PtclBdry(int         total_steps  =50000,             // total_steps
		      int         snapNum    =100,               // number of snapNum   
		      int         interval     =10,                // print interval
		      REAL        dimn         =0.05,              // dimension of particle-composed-boundary
		      REAL        rsize        =1.0,               // relative container size
		      const char* iniptclfile  ="dep_particle_end",// input file, initial particles
		      const char* Particlefile ="scl_particle",    // output file, resulted particles, including snapNum 
		      const char* contactfile  ="scl_contact",     // output file, resulted contacts, including snapNum
		      const char* progressfile ="scl_progress",    // output file, statistical info
		      const char* debugfile    ="scl_debug");      // output file, debug info
  
  
  // generate particles in space for particle boundaries
  void generate_p(Gradation& grad,
		  const char* str,
		  int freetype,
		  REAL rsize,
		  REAL ht);
  
 
  void deGravitation(int   total_steps,  
		     int   snapNum,
		     int   interval,
		     bool  toRebuild,
		     const char* iniptclfile,   
		     const char* Particlefile, 
		     const char* contactfile,
		     const char* progressfile, 
		     const char* debugfile);
  
  // actual deposit function for particle boundaries
  void deposit_p(int         total_steps  =50000,             // total_steps
		 int         snapNum    =100,               // number of snapNum   
		 int         interval     =10,                // print interval 
		 REAL dimn   =0.05,                           // dimension of particle-composed-boundary
		 REAL rsize  =1.0,                            // relative container size
		 const char* iniptclfile  ="flo_particle_end",// input file, initial particles
		 const char* Particlefile ="dep_particle",    // output file, resulted particles, including snapNum 
		 const char* contactfile  ="dep_contact",     // output file, resulted contacts, including snapNum
		 const char* progressfile ="dep_progress",    // output file, statistical info
		 const char* debugfile    ="dep_debug");      // output file, debug info
  
  //squeeze paticles inside a container by moving the boundaries
  void squeeze(int         total_steps  =20000,               // total_steps
	       int         init_steps   =5000,                // initial_steps to reach equilibrium
	       int         snapNum    =100,                 // number of snapNum   
	       int         interval     =10,                  // print interval 
	       int         flag         =-1,                  // -1 squeeze; +1 loosen
	       const char* iniptclfile  ="flo_particle_end",  // input file, initial particles
	       const char* inibdryfile  ="dep_boundary_ini",  // input file, initial boundaries
	       const char* Particlefile ="dep_particle",      // output file, resulted particles, including snapNum 
	       const char* boundaryfile ="dep_boundary",      // output file, resulted boundaries
	       const char* contactfile  ="dep_contact",       // output file, resulted contacts, including snapNum
	       const char* progressfile ="dep_progress",      // output file, statistical info
	       const char* debugfile    ="dep_debug");        // output file, debug info
  
  void deposit_repose(int   interval,
		      const char* inibdryfile,
		      const char* Particlefile, 
		      const char* contactfile,
		      const char* progressfile, 
		      const char* debugfile);
  
  void angleOfRepose(int   interval,
		     const char* inibdryfile,
		     const char* Particlefile, 
		     const char* contactfile,
		     const char* progressfile, 
		     const char* debugfile);
  
  REAL getPtclMinX(const std::vector<Particle*> &particleVec) const;
  REAL getPtclMaxX(const std::vector<Particle*> &particleVec) const;
  REAL getPtclMinY(const std::vector<Particle*> &particleVec) const;
  REAL getPtclMaxY(const std::vector<Particle*> &particleVec) const;
  REAL getPtclMinZ(const std::vector<Particle*> &particleVec) const;
  REAL getPtclMaxZ(const std::vector<Particle*> &particleVec) const;
  
  void collapse(int   total_steps,  
		int   snapNum,
		int   interval,
		const char* iniptclfile,
		const char* initboundary,
		const char* Particlefile,
		const char* contactfile,
		const char* progressfile,
		const char* debugfile);
  
  void createMemParticle(REAL rRadius,
			 bool toRebuild,
			 const char* Particlefile,
			 const char* allParticle);
  
  void iso_MemBdry(int   total_steps,  
		   int   snapNum, 
		   int   interval,
		   REAL  sigma3,
		   REAL  rRadius,
		   bool  toRebuild,
		   const char* iniptclfile, 
		   const char* Particlefile,
		   const char* contactfile, 
		   const char* progressfile,
		   const char* debugfile);
  
  void TrimPtclBdryByHeight(REAL height,
			    const char* iniptclfile,
			    const char* Particlefile);
  
  void applyParticleBoundary(int          total_steps  =100000,
			     int          snapNum    =100,
			     int          nterval      =10,
			     REAL         sigma        =1.0e+4,
			     const char*  iniptclfile  ="cre_particle",
			     const char*  inibdryfile  ="cre_bounary",
			     const char*  Particlefile ="iso_particle",
			     const char*  boundaryfile ="iso_boundary",
			     const char*  contactfile  ="iso_contact",
			     const char*  progressfile ="iso_progress",
			     const char*  balancedfile ="iso_balanced",
			     const char*  debugfile    ="iso_debug");
  
  // Isotropically compress floating particles to a specific confining pressure, which is usually a low
  // value in order to create an intial status. Force boundaries are used. This process may be not 
  // physically true.
  void isotropic(int          total_steps  =100000,
		 int          snapNum    =100,
		 int          interval     =10,
		 REAL         sigma        =1.0e+4,
		 const char*  iniptclfile  ="flo_particle_end",
		 const char*  inibdryfile  ="iso_inbdry",
		 const char*  Particlefile ="iso_particle",
		 const char*  boundaryfile ="iso_boundary",
		 const char*  contactfile  ="iso_contact",
		 const char*  progressfile ="iso_progress",
		 const char*  balancedfile ="iso_balanced",
		 const char*  debugfile    ="iso_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // increases confining pressure step by step to sigma_b, making it possible to find equilibrium 
  // state where particle pressure equals confining pressure. Force boundaries are used.
  void isotropic(int          total_steps   =100000,
		 int          snapNum     =100,
		 int          interval      =10, 
		 REAL  sigma_a       =1.0e+4,
		 REAL  sigma_b       =1.0e+5,	
		 int    sigma_division      =100,	  
		 const char*  iniptclfile   ="iso_particle_10k",
		 const char*  inibdryfile   ="iso_boundary_10k",
		 const char*  Particlefile  ="iso_particle", 
		 const char*  boundaryfile  ="iso_boundary", 
		 const char*  contactfile   ="iso_contact",
		 const char*  progressfile  ="iso_progress",
		 const char*  balancedfile  ="iso_balanced", 
		 const char*  debugfile     ="iso_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // follows an unloading-reloading stress path. Force boundaries are used.
  void isotropic(int          total_steps,
		 int          snapNum,
		 int          interval,
		 int          sigma_points,			  
		 REAL  sigma_values[],
		 int          sigma_division=100,
		 const char*  iniptclfile   ="iso_particle_10k",
		 const char*  inibdryfile   ="iso_boundary_10k",
		 const char*  Particlefile  ="iso_particle", 
		 const char*  boundaryfile  ="iso_boundary", 
		 const char*  contactfile   ="iso_contact",
		 const char*  progressfile  ="iso_progress",
		 const char*  balancedfile  ="iso_balanced", 
		 const char*  debugfile     ="iso_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_3. This function
  // increases confining pressure step by step to sigma_1, thus making it possible to find out
  // balanced status where top & bottom particle pressure equals major principle stress. 
  // Side boundaries are fixed, top and bottom plates are force-controlled.
  void odometer(int          total_steps    =100000,
		int          snapNum      =100,
		int          interval       =10,
		REAL  sigma_3        =1.0e+4,
		REAL  sigma_1        =1.0e+5,
		int          sigma_division =100,		  
		const char*  iniptclfile    ="iso_particle_10k",
		const char*  inibdryfile    ="iso_boundary_10k",
		const char*  Particlefile   ="odo_particle", 
		const char*  boundaryfile   ="odo_boundary", 
		const char*  contactfile    ="odo_contact",
		const char*  progressfile   ="odo_progress",
		const char*  balancedfile   ="odo_balanced", 
		const char*  debugfile      ="odo_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_3. This function
  // increases confining pressure step by step to sigma_1, thus making it possible to find out
  // balanced status where top & bottom particle pressure equals major principle stress. 
  // Side boundaries are fixed, top and bottom plates are force-controlled. Unloading is applied.
  void odometer(int          total_steps,
		int          snapNum,
		int          interval,
		int          sigma_points,			  
		REAL  sigma_values[],
		int          sigma_division=100,		  
		const char*  iniptclfile   ="iso_particle_10k",
		const char*  inibdryfile   ="iso_boundary_10k",
		const char*  Particlefile  ="odo_particle", 
		const char*  boundaryfile  ="odo_boundary", 
		const char*  contactfile   ="odo_contact",
		const char*  progressfile  ="odo_progress",
		const char*  balancedfile  ="odo_balanced", 
		const char*  debugfile     ="odo_debug");
  
  // The confining pressure is 500kPa. This function initializes triaxial compression test.
  void triaxialPtclBdryIni(int          total_steps  =10000,
			   int          snapNum    =100,
			   int          interval     =10,
			   REAL         sigma        =5.0e+5,
			   const char*  iniptclfile  ="ini_particle_ini",
			   const char*  inibdryfile  ="ini_boundary_ini",
			   const char*  Particlefile ="ini_particle", 
			   const char*  boundaryfile ="ini_boundary", 
			   const char*  contactfile  ="ini_contact",
			   const char*  progressfile ="ini_progress",
			   const char*  debugfile    ="ini_debug");
  
  // The confining pressure is 500kPa. This function performs triaxial compression test.
  // Displacement boundaries are used in axial direction.
  void triaxialPtclBdry(int          total_steps  =100000,
			int          snapNum    =100,
			int          interval     =10,
			const char*  iniptclfile  ="iso_particle_100k",
			const char*  inibdryfile  ="iso_boundary_100k",
			const char*  Particlefile ="tri_particle", 
			const char*  boundaryfile ="tri_boundary", 
			const char*  contactfile  ="tri_contact",
			const char*  progressfile ="tri_progress",
			const char*  balancedfile ="tri_balanced", 
			const char*  debugfile    ="tri_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // performs triaxial compression test. Displacement boundaries are used in axial direction.
  void triaxial(int          total_steps  =100000,
		int          snapNum    =100,
		int          interval     =10,
		REAL  sigma_a      =1.0e+5,
		const char*  iniptclfile  ="iso_particle_100k",
		const char*  inibdryfile  ="iso_boundary_100k",
		const char*  Particlefile ="tri_particle", 
		const char*  boundaryfile ="tri_boundary", 
		const char*  contactfile  ="tri_contact",
		const char*  progressfile ="tri_progress",
		const char*  balancedfile ="tri_balanced", 
		const char*  debugfile    ="tri_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // performs triaxial compression test with unloading. Displacement boundaries are used in 
  // axial direction.
  void triaxial(int          total_steps  =200000,
		int          unload_step  =100000,
		int          snapNum    =100,
		int          interval     =10,
		REAL  sigma_a      =3.0e+5,
		const char*  iniptclfile  ="iso_particle_300k",
		const char*  inibdryfile  ="iso_boundary_300k",
		const char*  Particlefile ="tri_particle", 
		const char*  boundaryfile ="tri_boundary", 
		const char*  contactfile  ="tri_contact",
		const char*  progressfile ="tri_progress",
		const char*  balancedfile ="tri_balanced", 
		const char*  debugfile    ="tri_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // A rectangular pile is then drived into the particles using displacement control.
  void rectPile_Disp(int          total_steps  =50000,
		     int          snapNum    =100,
		     int          interval     =10,
		     const char*  iniptclfile  ="pile_particle_ini",
		     const char*  inibdryfile  ="pile_boundary_ini",
		     const char*  Particlefile ="pile_particle", 
		     const char*  boundaryfile ="pile_boundary", 
		     const char*  contactfile  ="pile_contact",
		     const char*  progressfile ="pile_progress",
		     const char*  debugfile    ="pile_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // An ellipsoidal pile is then drived into the particles using displacement control.
  void ellipPile_Disp(int         total_steps  =50000,  
		      int         snapNum    =100, 
		      int          interval     =10,
		      REAL dimn         =0.05,
		      REAL rsize        =1.0,
		      const char* iniptclfile  ="pile_particle_ini",
		      const char* Particlefile ="pile_particle", 
		      const char* contactfile  ="pile_contact",  
		      const char* progressfile ="pile_progress",
		      const char* debugfile    ="pile_debug");
  
  // The specimen has been deposited with gravitation within rigid boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
  void ellipPile_Impact(int         total_steps  =50000,  
			int         snapNum    =100, 
			int         interval     =10,
			REAL dimn         =0.05,
			const char* iniptclfile  ="ipt_particle_ini",
			const char* inibdryfile  ="dep_boundary_ini",
			const char* Particlefile ="ipt_particle", 
			const char* contactfile  ="ipt_contact",  
			const char* progressfile ="ipt_progress",
			const char* debugfile    ="ipt_debug");
  
  // The specimen has been deposited with gravitation within particle boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
  void ellipPile_Impact_p(int         total_steps  =50000,  
			  int         snapNum    =100, 
			  int         interval     =10,
			  REAL dimn         =0.05,
			  const char* iniptclfile  ="ipt_particle_ini",
			  const char* Particlefile ="ipt_particle", 
			  const char* contactfile  ="ipt_contact",  
			  const char* progressfile ="ipt_progress",
			  const char* debugfile    ="ipt_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // An ellipsoidal pile is then drived into the particles using force control.
  void ellipPile_Force(int         total_steps  =50000,  
		       int         snapNum    =100, 
		       int         interval     =10,
		       REAL dimn         =0.05,
		       REAL force        =1.0e+4,
		       int   division           =100,
		       const char* iniptclfile  ="pile_particle_ini",
		       const char* Particlefile ="pile_particle", 
		       const char* contactfile  ="pile_contact",  
		       const char* progressfile ="pile_progress",
		       const char* balancedfile ="pile_balanced",
		       const char* debugfile    ="pile_debug");

  void truetriaxial(int          total_steps   =1000000,
		    int          snapNum     =100,
		    int          interval      =10,
		    REAL  sigma_a       =1.0e+4,
		    REAL  sigma_w       =1.0e+5,
		    REAL  sigma_l       =1.0e+5,	
		    REAL  sigma_h       =1.0e+5,	
		    int          sigma_division=100,			  
		    const char*  iniptclfile   ="iso_particle_10k",
		    const char*  inibdryfile   ="iso_boundary_10k",
		    const char*  Particlefile  ="tru_particle", 
		    const char*  boundaryfile  ="tru_boundary", 
		    const char*  contactfile   ="tru_contact",
		    const char*  progressfile  ="tru_progress",
		    const char*  balancedfile  ="tru_balanced", 
		    const char*  debugfile     ="tru_debug");

  public:
    void findParticleInRectangle(const Rectangle &container,
				 const std::vector<Particle*> &allParticle,
				 std::vector<Particle*> &foundParticle);
    void mpiTest0();
    
    void mpiTest1() {
      int tag = 0;
      if (boostWorld.rank() == 0) {
	readParticle("test_particle_ini");
	printParticle("test_particle_0");
	boostWorld.send(1, tag, allParticleVec);
	boostWorld.send(1, tag, gradation);
      }
      else if (boostWorld.rank() == 1 ) {
	boostWorld.recv(0, tag, allParticleVec);
	boostWorld.recv(0, tag, gradation);
	printParticle("test_particle_1");
      }
    }
    
    void mpiTest2() {
      if (boostWorld.rank() == 0) {
	readParticle("test_particle_ini");
	printParticle("test_particle_0");
      }
      broadcast(boostWorld, allParticleVec, 0);
      broadcast(boostWorld, gradation, 0);
      
      if (boostWorld.rank() == 1 ) {
	printParticle("test_particle_1");
      }
    }
    
  };
  
} // namespace dem

#endif
