#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "realtypes.h"
#include "vec.h"
#include "gradation.h"
#include "particle.h"
#include "contact.h"
#include "boundary.h"
#include "rectangle.h"
#include "cylinder.h"
#include "spring.h"
#include "matrix.h"
#include "cell.h"
#include "edge.h"
#include "microrve.h"
#include "omegabndry.h"
#include <map>
#include <vector>
#include <fstream>

namespace dem {

class assembly{
  typedef contact<particle>  CONTACT;
  typedef rgd_bdry<particle> RGDBDRY;
  typedef flb_bdry<particle> FLBBDRY;
  
 public:
  assembly(){
    TotalNum = 0;
    trimHistoryNum = 0;
    Volume = 0;
    BulkDensity = 0;
    PossCntctNum = 0;
    ActualCntctNum = 0;
    RgdBdryNum = 0;
    FlbBdryNum = 0;
  }
  
  ~assembly(){
    std::vector<particle*>::iterator pt;
    std::vector<RGDBDRY*>::iterator  rt;
    std::vector<FLBBDRY*>::iterator  ft;
    std::vector<spring*>::iterator st;
    std::vector<cell*>::iterator ct;

    std::vector<microrve*>::iterator mt;
    // it is important to release memory pointed to by pointers in the container,
    // otherwise memory leaking occurs
    for(pt = ParticleVec.begin(); pt != ParticleVec.end(); ++pt)
      delete (*pt);
    for(rt = RBVec.begin(); rt != RBVec.end(); ++rt)
      delete (*rt);
    for(ft = FBVec.begin(); ft!=FBVec.end(); ++ft)
      delete (*ft);
    for(st = SpringVec.begin(); st != SpringVec.end(); ++st)
      delete (*st);
    for(ct = cellVec.begin(); ct != cellVec.end(); ++ct)
      delete (*ct);
    for(ct = cellVec_init.begin(); ct != cellVec_init.end(); ++ct)
      delete (*ct);

    for(mt = RVEVec.begin(); mt != RVEVec.end(); ++mt)
      delete (*mt);
    
    // in case of consecutive simulations
    ParticleVec.clear();
    RBVec.clear();
    FBVec.clear();
    SpringVec.clear();
    cellVec.clear();
    cellVec_init.clear();
    RVEVec.clear();
  }
 
  void setContainer(rectangle &cont) {container = cont;} 
  void setGradation(gradation &grad) {gradInfo = grad;}
  void setCavity(rectangle &cavi) {cavity = cavi;}
  void readSample(const char* str);           // create a sample with particles from an existing file
  void readRigidBoundary(std::ifstream &ifs); // create rigid boundaries from an existing file
  void readFlexiBoundary(std::ifstream &ifs); // create flxible boundaries from an existing file
  void readBoundary(const char* str);         // create either rigid or flexible boundaries from an existing file
  void trimCavity(bool toRebuild, const char* particlefile, const char* cavparticle);
  void readCavityBoundary(const char* boundaryfile);
  void buildCavityBoundary(int existMaxId, const char* boundaryfile);
  void findContact();                           // detect and resolve contact between particles
  void findParticleOnBoundary();                // find particles on boundaries
  void findParticleOnCavity();                  // find particle on cavity boundaries
  void findParticleOnLine();                    // find particles on lines
  void createFlbNet();
  
  void clearForce();                            // clear forces and moments for all particles
  void flexiBoundaryForceZero();
  void initFBForce();
  void internalForce(REAL& avgnm, REAL& avgsh); // calculate inter-particle forces
  void springForce();
  void rigidBoundaryForce();                    // calcualte forces between rigid boundaries and particles
  void rigidBoundaryForce(REAL penetr[],int cntnum[]);
  void flexiBoundaryForce();
  void cavityBoundaryForce();
  void updateParticle();                        // update motion of particles
  
  REAL ellipPileForce();                        // for force pile only
  void ellipPileUpdate();                       // for force pile only
  
  vec  ellipPileDimn();
  REAL ellipPileTipZ();
  REAL ellipPilePeneVol();
  
  // if bn[i]=2, the 2nd rigid boundary should be updated according to rbctl[i],
  // totally num rigid boundaries must be updated
  void updateRB(int bn[], UPDATECTL rbctl[], int num);
  void updateRB6();
  void updateRectPile();
  
  // if bn[i]=2, the 2nd flxible boundary should be updated according to fbctl[bn[i]*2-2] and fbctl[bn[i]*2-1], 
  // totally num flxible boundaries must be updated, the size of fbctl is 2 times large as size of bn
  void updateFB(int bn[], UPDATECTL fbctl[], int num);
  
  REAL getDensity() const; 
  int  getPossCntctNum() const {return  PossCntctNum;};
  int  getActualCntctNum() const {return ActualCntctNum;}
  REAL getAveragePenetration() const;
  REAL getVibraTimeStep() const;
  REAL getImpactTimeStep() const;
  REAL getAverageVelocity() const;
  REAL getAverageForce() const;
  REAL getAverageOmga() const;
  REAL getAverageMoment() const;
  REAL getParticleVolume() const;
  vec  getTopFreeParticlePosition() const;
  REAL getTransEnergy() const;
  REAL getRotatEnergy() const;
  REAL getKinetEnergy() const;
  REAL getPotenEnergy(REAL ref) const;
  
  matrix getGranularStress() const;	// calculate the granular stress, written on Feb 13, 2013, return a 3x3 matrix
  matrix getFabric() const;	// calculate the fabric tensor, written on March 8, 2013, return a 3x3 matrix. unit is m^2
  // for granular strain calculation, written on Feb 18, 2013
  void callQhull() const;	// based on input_for_Qhull that is from createInputForQhull() to tessellate and output "tess_info"
  void resetStartPosition();	// reset the initial position of each particle. For granular strain March 18, 2013
  void setThreshold(REAL threshold) { strainThreshold = threshold;}
  matrix getum(int ID, int type) const;	// calculate the displacement vector of particle ID, retur a 3x1 matrix
  matrix getb(cell, int ID) const;	// calculate the b vector of a cell, return a 3x1 matrix
  matrix getdmn(edge, std::vector<cell>) const;	// calculate the complementary area vector of edge
  matrix getGranularStrain();	// calculate granular strain, not just granular e
  void createInputForQhull() const;	// create input file for Qhull. March 18, 2013
  void readTesse(const char* str);	// create tessellation information from Qhull output files, finally it will not be used siince we need to call Qhull within the ellip3d code
  // for finite granular strain calculation, March 26, 2013
  void readTesse_finite(const char* str);	// create cell information into std::vector<cell> cellVec from Qhull out put files
  void setNumberingOrder();	// ensure that the numbering of each tetrahedron in std::vector<cell> cellVec is couter-clockwise
  matrix getBigB(const cell&) const;
  matrix getdudx(const cell&) const;	// calculate the matrix dudx, which is for F=dudx+I
  matrix getdudx_curr(const cell&) const;	// calculate the matrix dudx with respect to current configuration
  REAL getCellVolume(const cell&) const;	// calculate the volume of cell
  matrix getFiniteStrain();	// calculate finite granular strain using FEM-similar way
  matrix getHigherStrain() const;	// calculate 1/(2V)*sum(dudx'*dudx*vL), for mixed finite strain
  matrix getEulerianStrain();	// calculate Eulerian finite strain
  matrix getEuler_HOT()	const;
  matrix getdvdx_curr(const cell&) const;	// calculate spatial velocity gradient tensor for each tet, June 24, 2013
  matrix getAverage_dvdx() const;	// calculate average spatial velocity gradient tensor, June 24, 2013

  vec  getNormalForce(int bdry) const;       // get normal force acting on the bdry_th rigid boundary
  vec  getShearForce(int bdry) const;        // get shear force acting on the bdry_th rigid boundary
  REAL getAvgNormal(int bdry) const;
  vec  getApt(int bdry) const;               // get a point on bdry_th rigid boundary
  vec  getDirc(int bdry) const;              // get the dirc of bdry_th rigid boundry
  REAL getArea(int bdry) const;
  REAL getAverageRigidPressure() const;
  void setArea(int bdry,REAL a);             // set the area of the bdry-th rigid boundary be a
  void setTrimHistoryNum(int n) {trimHistoryNum = n;}
  void printParticle(const char* str) const; // print particles info into a disk file
  void printMemParticle(const char* str) const; // print membrane particles info into a disk file
  void plotSpring(const char *str) const;    // print springs in Tecplot format
  void plotBoundary(const char *str) const;
  void plotCavity(const char *str) const;
  void checkMembrane(std::vector<REAL> &vx ) const;
  void printContact(const char* str) const;  // print contacts information
  void printBoundary(const char* str) const; // print rigid boundaries info to a disk file
  void printRectPile(const char* str);       // append rectangular pile info into a disk file
  void printCavityBoundary(const char* str) const; // print cavity boundaries
  void printCavityParticle(int total, const char* str) const;

  void expandCavityParticles(bool toRebuild,
			     REAL percent,
			     const char* cavityptclfile,
			     const char* particlefile,
			     const char* newptclfile);

  // create a specimen by depositing particles into rigid boundaries
  void deposit_RgdBdry(int   freetype,
		       int   total_steps,  
		       int   snapshots,
		       int   interval,
		       REAL  rFloHeight,
//		       REAL sampleVolumeRatio,   // Modified by Boning Zhang
		       const char* iniptclfile,   
		       const char* inibdryfile,
		       const char* particlefile, 
		       const char* contactfile,
		       const char* progressfile, 
		       const char* creparticle,
		       const char* creboundary,
		       const char* debugfile);
  
  // continue to deposit after a cavity is created inside the particle assemblage
  void depositAfterCavity(int   total_steps,  
			  int   snapshots,
			  int   interval,
			  const char* iniptclfile,   
			  const char* inibdryfile,
			  const char* inicavefile,
			  const char* particlefile, 
			  const char* contactfile,
			  const char* progressfile, 
			  const char* debugfile);

  // create a specimen by depositing particles into particle boundaries
  void deposit_PtclBdry(gradation& grad,
			int   freetype,
			REAL rsize,
			int   total_steps,  
			int   snapshots,
			int   interval,
			const char* iniptclfile,   
			const char* particlefile, 
			const char* contactfile,
			const char* progressfile, 
			const char* debugfile);
  
  // scale the assembly with particle boundaries from deposited state until it reaches steady state
  void scale_PtclBdry(int         total_steps  =50000,             // total_steps
		      int         snapshots    =100,               // number of snapshots   
		      int         interval     =10,                // print interval
		      REAL        dimn         =0.05,              // dimension of particle-composed-boundary
		      REAL        rsize        =1.0,               // relative container size
		      const char* iniptclfile  ="dep_particle_end",// input file, initial particles
		      const char* particlefile ="scl_particle",    // output file, resulted particles, including snapshots 
		      const char* contactfile  ="scl_contact",     // output file, resulted contacts, including snapshots
		      const char* progressfile ="scl_progress",    // output file, statistical info
		      const char* debugfile    ="scl_debug");      // output file, debug info
  
  // generate particles in space for rigid boundaries
  void generate(const char* str,
		int freetype,
		REAL ht);
//     void generate(const char* str,
// 		int freetype,
//  		REAL vr);         // Modified by Boning Zhang

  
  // generate particles in space for particle boundaries
  void generate_p(gradation& grad,
		  const char* str,
		  int freetype,
		  REAL rsize,
		  REAL ht);
  
  // actual deposit function for rigid boundaries
  void deposit(int         total_steps  =100000,              // total_steps
	       int         snapshots    =100,                 // number of snapshots   
	       int         interval     =10,                  // print interval 
	       const char* iniptclfile  ="flo_particle_end",  // input file, initial particles
	       const char* inibdryfile  ="dep_boundary_ini",  // input file, initial boundaries
	       const char* particlefile ="dep_particle",      // output file, resulted particles, including snapshots 
	       const char* contactfile  ="dep_contact",       // output file, resulted contacts, including snapshots
	       const char* progressfile ="dep_progress",      // output file, statistical info
	       const char* debugfile    ="dep_debug");        // output file, debug info
  
  void deGravitation(int   total_steps,  
		     int   snapshots,
		     int   interval,
		     bool  toRebuild,
		     const char* iniptclfile,   
		     const char* particlefile, 
		     const char* contactfile,
		     const char* progressfile, 
		     const char* debugfile);
  
  // actual deposit function for particle boundaries
  void deposit_p(int         total_steps  =50000,             // total_steps
		 int         snapshots    =100,               // number of snapshots   
		 int         interval     =10,                // print interval 
		 REAL dimn   =0.05,                           // dimension of particle-composed-boundary
		 REAL rsize  =1.0,                            // relative container size
		 const char* iniptclfile  ="flo_particle_end",// input file, initial particles
		 const char* particlefile ="dep_particle",    // output file, resulted particles, including snapshots 
		 const char* contactfile  ="dep_contact",     // output file, resulted contacts, including snapshots
		 const char* progressfile ="dep_progress",    // output file, statistical info
		 const char* debugfile    ="dep_debug");      // output file, debug info
  
  //squeeze paticles inside a container by moving the boundaries
  void squeeze(int         total_steps  =20000,               // total_steps
	       int         init_steps   =5000,                // initial_steps to reach equilibrium
	       int         snapshots    =100,                 // number of snapshots   
	       int         interval     =10,                  // print interval 
	       int         flag         =-1,                  // -1 squeeze; +1 loosen
	       const char* iniptclfile  ="flo_particle_end",  // input file, initial particles
	       const char* inibdryfile  ="dep_boundary_ini",  // input file, initial boundaries
	       const char* particlefile ="dep_particle",      // output file, resulted particles, including snapshots 
	       const char* boundaryfile ="dep_boundary",      // output file, resulted boundaries
	       const char* contactfile  ="dep_contact",       // output file, resulted contacts, including snapshots
	       const char* progressfile ="dep_progress",      // output file, statistical info
	       const char* debugfile    ="dep_debug");        // output file, debug info
  
  void deposit_repose(int   interval,
		      const char* inibdryfile,
		      const char* particlefile, 
		      const char* contactfile,
		      const char* progressfile, 
		      const char* debugfile);
  
  void angleOfRepose(int   interval,
		     const char* inibdryfile,
		     const char* particlefile, 
		     const char* contactfile,
		     const char* progressfile, 
		     const char* debugfile);
  
  REAL getMaxCenterHeight() const;
  
  void collapse(int   total_steps,  
		int   snapshots,
		int   interval,
		const char* iniptclfile,
		const char* initboundary,
		const char* particlefile,
		const char* contactfile,
		const char* progressfile,
		const char* debugfile);
  
  void buildBoundary(int bdrynum, const char* boundaryfile);
  
  void buildBoundary(const char* boundaryfile);
  
  void trim(bool  toRebuild,
	    const char* particlefile,
	    const char* trmparticle);
  
  void createMemParticle(REAL rRadius,
			 bool toRebuild,
			 const char* particlefile,
			 const char* allparticle);
  
  void iso_MemBdry(int   total_steps,  
		   int   snapshots, 
		   int   interval,
		   REAL  sigma3,
		   REAL  rRadius,
		   bool  toRebuild,
		   const char* iniptclfile, 
		   const char* particlefile,
		   const char* contactfile, 
		   const char* progressfile,
		   const char* debugfile);
  
  void TrimPtclBdryByHeight(REAL height,
			    const char* iniptclfile,
			    const char* particlefile);
  
  void applyParticleBoundary(int          total_steps  =100000,
			     int          snapshots    =100,
			     int          nterval      =10,
			     REAL         sigma        =1.0e+4,
			     const char*  iniptclfile  ="cre_particle",
			     const char*  inibdryfile  ="cre_bounary",
			     const char*  particlefile ="iso_particle",
			     const char*  boundaryfile ="iso_boundary",
			     const char*  contactfile  ="iso_contact",
			     const char*  progressfile ="iso_progress",
			     const char*  balancedfile ="iso_balanced",
			     const char*  debugfile    ="iso_debug");
  
  // Isotropically compress floating particles to a specific confining pressure, which is usually a low
  // value in order to create an intial status. Force boundaries are used. This process may be not 
  // physically true.
  void isotropic(int          total_steps  =100000,
		 int          snapshots    =100,
		 int          interval     =10,
		 REAL         sigma        =1.0e+4,
		 const char*  iniptclfile  ="flo_particle_end",
		 const char*  inibdryfile  ="iso_inbdry",
		 const char*  particlefile ="iso_particle",
		 const char*  boundaryfile ="iso_boundary",
		 const char*  contactfile  ="iso_contact",
		 const char*  progressfile ="iso_progress",
		 const char*  balancedfile ="iso_balanced",
		 const char*  debugfile    ="iso_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // increases confining pressure step by step to sigma_b, making it possible to find equilibrium 
  // state where particle pressure equals confining pressure. Force boundaries are used.
  void isotropic(int          total_steps   =100000,
		 int          snapshots     =100,
		 int          interval      =10, 
		 REAL  sigma_a       =1.0e+4,
		 REAL  sigma_b       =1.0e+5,	
		 int    sigma_division      =100,	  
		 const char*  iniptclfile   ="iso_particle_10k",
		 const char*  inibdryfile   ="iso_boundary_10k",
		 const char*  particlefile  ="iso_particle", 
		 const char*  boundaryfile  ="iso_boundary", 
		 const char*  contactfile   ="iso_contact",
		 const char*  progressfile  ="iso_progress",
		 const char*  balancedfile  ="iso_balanced", 
		 const char*  debugfile     ="iso_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // follows an unloading-reloading stress path. Force boundaries are used.
  void isotropic(int          total_steps,
		 int          snapshots,
		 int          interval,
		 int          sigma_points,			  
		 REAL  sigma_values[],
		 int          sigma_division=100,
		 const char*  iniptclfile   ="iso_particle_10k",
		 const char*  inibdryfile   ="iso_boundary_10k",
		 const char*  particlefile  ="iso_particle", 
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
		int          snapshots      =100,
		int          interval       =10,
		REAL  sigma_3        =1.0e+4,
		REAL  sigma_1        =1.0e+5,
		int          sigma_division =100,		  
		const char*  iniptclfile    ="iso_particle_10k",
		const char*  inibdryfile    ="iso_boundary_10k",
		const char*  particlefile   ="odo_particle", 
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
		int          snapshots,
		int          interval,
		int          sigma_points,			  
		REAL  sigma_values[],
		int          sigma_division=100,		  
		const char*  iniptclfile   ="iso_particle_10k",
		const char*  inibdryfile   ="iso_boundary_10k",
		const char*  particlefile  ="odo_particle", 
		const char*  boundaryfile  ="odo_boundary", 
		const char*  contactfile   ="odo_contact",
		const char*  progressfile  ="odo_progress",
		const char*  balancedfile  ="odo_balanced", 
		const char*  debugfile     ="odo_debug");
  
  // The confining pressure is 500kPa. This function initializes triaxial compression test.
  void triaxialPtclBdryIni(int          total_steps  =10000,
			   int          snapshots    =100,
			   int          interval     =10,
			   REAL         sigma        =5.0e+5,
			   const char*  iniptclfile  ="ini_particle_ini",
			   const char*  inibdryfile  ="ini_boundary_ini",
			   const char*  particlefile ="ini_particle", 
			   const char*  boundaryfile ="ini_boundary", 
			   const char*  contactfile  ="ini_contact",
			   const char*  progressfile ="ini_progress",
			   const char*  debugfile    ="ini_debug");
  
  // The confining pressure is 500kPa. This function performs triaxial compression test.
  // Displacement boundaries are used in axial direction.
  void triaxialPtclBdry(int          total_steps  =100000,
			int          snapshots    =100,
			int          interval     =10,
			const char*  iniptclfile  ="iso_particle_100k",
			const char*  inibdryfile  ="iso_boundary_100k",
			const char*  particlefile ="tri_particle", 
			const char*  boundaryfile ="tri_boundary", 
			const char*  contactfile  ="tri_contact",
			const char*  progressfile ="tri_progress",
			const char*  balancedfile ="tri_balanced", 
			const char*  debugfile    ="tri_debug");
  
  // The specimen has been isotropically compressed to confining pressure sigma_a. This function
  // performs triaxial compression test. Displacement boundaries are used in axial direction.
  void triaxial(int          total_steps  =100000,
		int          snapshots    =100,
		int          interval     =10,
		REAL  sigma_a      =1.0e+5,
		const char*  iniptclfile  ="iso_particle_100k",
		const char*  inibdryfile  ="iso_boundary_100k",
		const char*  particlefile ="tri_particle", 
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
		int          snapshots    =100,
		int          interval     =10,
		REAL  sigma_a      =3.0e+5,
		const char*  iniptclfile  ="iso_particle_300k",
		const char*  inibdryfile  ="iso_boundary_300k",
		const char*  particlefile ="tri_particle", 
		const char*  boundaryfile ="tri_boundary", 
		const char*  contactfile  ="tri_contact",
		const char*  progressfile ="tri_progress",
		const char*  balancedfile ="tri_balanced", 
		const char*  debugfile    ="tri_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // A rectangular pile is then drived into the particles using displacement control.
  void rectPile_Disp(int          total_steps  =50000,
		     int          snapshots    =100,
		     int          interval     =10,
		     const char*  iniptclfile  ="pile_particle_ini",
		     const char*  inibdryfile  ="pile_boundary_ini",
		     const char*  particlefile ="pile_particle", 
		     const char*  boundaryfile ="pile_boundary", 
		     const char*  contactfile  ="pile_contact",
		     const char*  progressfile ="pile_progress",
		     const char*  debugfile    ="pile_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // An ellipsoidal pile is then drived into the particles using displacement control.
  void ellipPile_Disp(int         total_steps  =50000,  
		      int         snapshots    =100, 
		      int          interval     =10,
		      REAL dimn         =0.05,
		      REAL rsize        =1.0,
		      const char* iniptclfile  ="pile_particle_ini",
		      const char* particlefile ="pile_particle", 
		      const char* contactfile  ="pile_contact",  
		      const char* progressfile ="pile_progress",
		      const char* debugfile    ="pile_debug");
  
  // The specimen has been deposited with gravitation within rigid boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
  void ellipPile_Impact(int         total_steps  =50000,  
			int         snapshots    =100, 
			int         interval     =10,
			REAL dimn         =0.05,
			const char* iniptclfile  ="ipt_particle_ini",
			const char* inibdryfile  ="dep_boundary_ini",
			const char* particlefile ="ipt_particle", 
			const char* contactfile  ="ipt_contact",  
			const char* progressfile ="ipt_progress",
			const char* debugfile    ="ipt_debug");
  
  // The specimen has been deposited with gravitation within particle boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
  void ellipPile_Impact_p(int         total_steps  =50000,  
			  int         snapshots    =100, 
			  int         interval     =10,
			  REAL dimn         =0.05,
			  const char* iniptclfile  ="ipt_particle_ini",
			  const char* particlefile ="ipt_particle", 
			  const char* contactfile  ="ipt_contact",  
			  const char* progressfile ="ipt_progress",
			  const char* debugfile    ="ipt_debug");
  
  // The specimen has been deposited with gravitation within boundaries composed of particles.
  // An ellipsoidal pile is then drived into the particles using force control.
  void ellipPile_Force(int         total_steps  =50000,  
		       int         snapshots    =100, 
		       int         interval     =10,
		       REAL dimn         =0.05,
		       REAL force        =1.0e+4,
		       int   division           =100,
		       const char* iniptclfile  ="pile_particle_ini",
		       const char* particlefile ="pile_particle", 
		       const char* contactfile  ="pile_contact",  
		       const char* progressfile ="pile_progress",
		       const char* balancedfile ="pile_balanced",
		       const char* debugfile    ="pile_debug");

  void truetriaxial(int          total_steps   =1000000,
		    int          snapshots     =100,
		    int          interval      =10,
		    REAL  sigma_a       =1.0e+4,
		    REAL  sigma_w       =1.0e+5,
		    REAL  sigma_l       =1.0e+5,	
		    REAL  sigma_h       =1.0e+5,	
		    int          sigma_division=100,			  
		    const char*  iniptclfile   ="iso_particle_10k",
		    const char*  inibdryfile   ="iso_boundary_10k",
		    const char*  particlefile  ="tru_particle", 
		    const char*  boundaryfile  ="tru_boundary", 
		    const char*  contactfile   ="tru_contact",
		    const char*  progressfile  ="tru_progress",
		    const char*  balancedfile  ="tru_balanced", 
		    const char*  debugfile     ="tru_debug");
  
  void unconfined(int          total_steps  =100000,
		  int          snapshots    =100,
		  int          interval     =10, 
		  const char*  iniptclfile  ="flo_particle_end",
		  const char*  inibdryfile  ="unc_inbdry",
		  const char*  particlefile ="unc_particle", 
		  const char*  contactfile  ="unc_contact",
		  const char*  progressfile ="unc_progress",
		  const char*  debugfile    ="unc_debug");
  
  void soft_tric(REAL _sigma3, REAL _b,
		 const char* iniptclfile ="isotropic",
		 const char* boundaryfile="isobdry",
		       const char* responsefile="sftc",
		 const char* resultfile  ="sftcompressed",
		 const char* trackfile   ="tracksft");
  
  void earthPressure(REAL pressure, bool IsPassive,
		     const char* iniptclfile ="isotropic",
		     const char* boundaryfile="isobdry",
		     const char* responsefile="ssvpp",
		     const char* resultfile  ="pp",
		     const char* trackfile   ="pptrack");

  void shallowFoundation(const char* iniptclfile ="isotropic",
			 const char* boundaryfile="isobdry",
			 const char* responsefile="ssvsf",
			 const char* resultfile  ="shallowcompressed",
			 const char* trackfile   ="shallowtrack");
  
  void simpleShear(REAL normal_pressure,REAL _b,
		   const char* iniptclfile ="isotropic",
		   const char* boundaryfile="isobdry",
		   const char* responsefile="simpleshear",
		   const char* resultfile  ="simplesheared",
		   const char* trackfile   ="tracksimple");
  
  void dircShear(REAL rate, REAL roterate, REAL stress,
		 const char* iniptclfile ="deposit",
		 const char* boundaryfile="rgdcube.data",
		 const char* responsefile="ssvds",
		 const char* resultfile  ="dircshear",
		 const char* trackfile   ="trackds");

// calculate granular stress for assembly from input files, such as the assembly after deposition.
// writtend on Feb 16, 2013
  void calculateStress(REAL   relativeHeight,	// ratio of height of each subdomain to boundary height (assume they are equal)
		       	       REAL   relativeWidth,	// ratio of width of each subdomain to boundary width
			       REAL   relativeLength,	// ratio of length of each subdomain to boundary length
			       std::vector<REAL> &relBaseH,	// ratio of the base height of each subdomain above the assembly boundary base
			       const char* inptclfile,	// particle file, can use trm_particle_end for deposited assembly
			       const char* inbdryfile,
			       const char* resultfile);

  void trimByHeight(REAL Height,
			const char* inptclefile,	// initial particle file
			const char* inbdryfile,		// initial boundary file
			const char* resultptcle,
			const char* resultbdry);	// result particle file
  
  // this will just calculate the stress at the middle height of the assembly, and the will calculate formula stress based on the relBaseH vector. March 4, 2013
  void calculateMiddleStress(int subNum,
			       const char* inptclfile,	// particle file, can use trm_particle_end for deposited assembly
			       const char* inbdryfile,
			       const char* resultfile);
// calculate the fabric tensor of an assembly based on input particle file, and then output each element of the fabric tensor and its three invariants. March 8, 2013
  void calculateFabric(const char* inptclfile,
		       const REAL boxL,	// length, width and height of box, used for porosity
                       const REAL boxW,	// unit is mm
		       const REAL boxH,
		       const char* resultfile);

// reset coordinates of particles which are contacting with boundaries. And then call Qhull to tessellate. This function is written for Nonlinear FEM term project -- a 3D tet poromechanics Matlab code.
  void resetCoord_FEM(const char* inptclfile,
		      const char* inbdryfile);

// rotate the boundary walls along different axis. Test finite granular strain. May 21, 2013
  void rotateXYZ(REAL angle,	// angle that you want to roate
	       vec axis_flag,	// represent axis you want to roate: 0 means not to this axis, 1 means along this axis. If axis_num is 1 1 1, then it means rotate along x,y,z
	       int          snapshots    =100,
	       int          interval     =10,
	       const char*  iniptclfile  ="flo_particle_end",
	       const char*  inibdryfile  ="iso_inbdry",
	       const char*  particlefile ="rotate_particle",
	       const char*  boundaryfile ="rotate_boundary",
	       const char*  contactfile  ="rotate_contact",
	       const char*  progressfile ="rotate_progress",
	       const char*  balancedfile ="rotate_balanced",
	       const char*  debugfile    ="rotate_debug");


// apply the same rotation tensor to each particle directly such that we can make sure that the rotations of all particles are the same pure rotation. Test finite granular strain. June 10, 2013
  void rotateXYZ_dir(REAL angle,	// angle that you want to roate
	       vec axis_flag,	// represent axis you want to roate: 0 means not to this axis, 1 means along this axis. If axis_num is 1 1 1, then it means rotate along x,y,z
	       int          snapshots    =100,
	       int          interval     =10,
	       const char*  iniptclfile  ="flo_particle_end",
	       const char*  inibdryfile  ="iso_inbdry",
	       const char*  particlefile ="rotate_particle",
	       const char*  boundaryfile ="rotate_boundary",
	       const char*  contactfile  ="rotate_contact",
	       const char*  progressfile ="rotate_progress",
	       const char*  balancedfile ="rotate_balanced",
	       const char*  debugfile    ="rotate_debug");



  // this is actually a deposit process, but since we cannot calculate granular stress and strain during normal deposition, however in cavity expansion we need to calculate granular stress and strain. Then we write a new function cavityExpand to do cavity expansion and in this function, we will calculate granular stress and strain. June 17, 2013
  void cavityExpand(int         total_steps  =100000,              // total_steps
	       int         snapshots    =100,                 // number of snapshots   
	       int         interval     =10,                  // print interval 
	       const char* iniptclfile  ="flo_particle_end",  // input file, initial particles
	       const char* inibdryfile  ="dep_boundary_ini",  // input file, initial boundaries
	       const char* particlefile ="exp_particle",      // output file, resulted particles, including snapshots 
	       const char* contactfile  ="exp_contact",       // output file, resulted contacts, including snapshots
	       const char* progressfile ="exp_progress",      // output file, statistical info
	       const char* debugfile    ="exp_debug");        // output file, debug info

  // this is used to calculate volume of the assembly based on particle input file, July 14, 2013
  void calculateVolume(const char* iniptclefile);


  // return granular stress in the domain, Jan 17, 2014
  matrix getGranularStress(REAL xmin, REAL xmax,	// dimensions of the domain in which 
			   REAL ymin, REAL ymax,	// we want to caclulate granular stress
			   REAL zmin, REAL zmax );



  // micromophic granular stress, Jan 7, 2014

  // create RVEVec and bndryRVEVec, 
  // Omega is the interesting domain in which we want to calculate micromophic granular stress
  void generateRVE(int nx, int ny, int nz,	// these are the numbers of RVEs in the three directions
		   REAL xOmega_min, REAL xOmega_max, 
		   REAL yOmega_min, REAL yOmega_max, 
		   REAL zOmega_min, REAL zOmega_max);								 

  void calcRVEUnSymmetricGranularStress();	// calculate the symmetric Granular stress, March 6 2014
  void calcRVErho_xi(bool flag);		// calculate xi_arm_nplus and xi_arm_n for each RVE
  void calcRVERelativePosition(bool flag);
  void calcRVEa_xi(REAL dt);	// calculate the acceleration of xi_arm
  // calculate the acceleration of the center of Omega, to check the balance of linear momentum, March 18, 2014
  void calca_Omega_center(REAL dt){a_Omega_center = (coord_Omega_center-2*coord_Omega_center_n+coord_Omega_center_nminus)/(dt*dt);}	

  void calcV_Omega();
  void initRVExi_n();		// set xi_arm_n = xi_arm_nplus, at initial
  void updateRVExi();		// set xi_arm_nminus = xi_arm_n & xi_arm_n = xi_arm_nplus for each RVE
  // set coord_Omega_center_n = coord_Omega_center, at initial, to check the balance of linear momentum, March 18, 2014
  void initOmega_center_n(){coord_Omega_center_n = coord_Omega_center;}	
  // set coord_Omega_center_nminus = coord_Omega_center_n & coord_Omega_center_n = coord_Omega_center, to check the balance of linear momentum, March 18, 2014
  void updateOmega_center(){coord_Omega_center_nminus = coord_Omega_center_n; coord_Omega_center_n = coord_Omega_center;}	
  void updateOmega_z(int , REAL, REAL);		// changing Omega and RVEs, update coord_Omega_center, min, max and v in RVE
						// update the boundary face area, 
						// z means Omega is changing in z direction for cavity expansion

  // March 17, 2014. update Omega and RVEs in x, y and z direction, for isotropic and triaxial compression
  void updateOmega_xyz(int nx, int ny, int nz, 
		       REAL xmin, REAL xmax,
		       REAL ymin, REAL ymax,
		       REAL zmin, REAL zmax  );	// changing Omega and RVEs, update coord_Omega_center, min, max and v in RVE
						// update the boundary face area, 
						// z means Omega is changing in z direction for cavity expansion

//  void calcXi_arm();  
  matrix getMicroS_mk();	// calculate micromophic symmetric granular stress based on grain stress in RVEs
  matrix getMicroM_lkm(int l);  // where l is 1, 2 or 3; l=1: microM_lkm as l=1
  matrix getMicroSigma(matrix &, matrix &); 	// the microSigma_acceOff, microSigma_bothOff
  matrix getLinearMomentumResidual();		// get the residual of linear momentum balance equation, March 18, 2014
  matrix calculateMicroDeformation();		// calculate micro-deformation tensor





 private:
  // particles property
  int  TotalNum;                      // total number of particles
  int  trimHistoryNum;                // record history maximum numbering before trimming
  int  PossCntctNum;                  // possible contact number based on spherical distances
  int  ActualCntctNum;                // actual contact number based on solution of 6th order equations
  std::vector<particle*> ParticleVec; // a vector of pointers, each pointing to a particle
  std::vector<CONTACT>   ContactVec;  // a vector of contacts
  std::vector<cnttgt>    CntTgtVec;   // a vector to store tangential contact force and displacement
  std::vector< std::vector< std::vector<particle*>  >  >   MemBoundary;  // a vector to store six boundaries
  std::vector<spring*>   SpringVec;   // a vector to store springs;
  gradation              gradInfo;    // particles gradation
  
  // container property
  rectangle container;
  rectangle cavity;
  REAL Volume;      // volume of the specimen
  REAL initVolume;  // initial volume of the specimen
  REAL BulkDensity; // bulk density of specimen
  
  REAL strainThreshold;	// for granular strain calculation, used to control tessellate.
  // boundary property
  int  BdryType;               // 0 - rigid boundaries; 1 - flxible boundaries
  int  RgdBdryNum;             // rigid boundary number
  int  FlbBdryNum;             // flxible boundary number
  std::vector<RGDBDRY*> RBVec; // a vector of pointers, each pointing to a rigid boundary surface
  std::vector<RGDBDRY*> CavityRBVec; // a vector of pointers, each pointing to a rigid cavity surface
  std::vector<FLBBDRY*> FBVec; // a vector of pointers, each pointing to a flexible boundary surface
  std::map<int,std::vector<boundarytgt>  > BdryTgtMap; // a map to store particle-boundary contacts' tangential info
  std::map<edge, std::vector<cell> > edge_map;	// used to store the edge set from tessellation, for granular strain calculation
  std::vector<cell*> cellVec;	// used to store cell info, for finite granular strain calculation, March 26, 2013
  std::vector<cell*> cellVec_init;	// the initial cell vector, for lagrangian finite strain

  matrix average_dudx_Bagi;	// used to test quadratic terms (dudx)^T*dudx, April 22, 2013
  matrix average_dudx_Lagrangian;
  matrix average_dudx_Eulerian;



  // micromohpic granular stress, Jan 7, 2014
  std::vector<microrve*> RVEVec;		// vector of all the RVEs
  std::vector<omegabndry> bndryFaceVec;	// vector of boundary faces of Omega
  REAL V_Omega;				// the volume of Omega
  matrix coord_Omega_center;		// center of Omega, row vector
  matrix coord_Omega_center_n;		// center of Omega, at previous step, i.e. step n
  matrix coord_Omega_center_nminus;	// center of Omega, at step n-1
  matrix a_Omega_center;		// acceleration of the center of Omega

//  matrix microS_mk;
//  matrix microM_1km;   			// microM_lkm as l=1
//  matrix microM_2km;			// microM_lkm as l=2
//  matrix microM_3km;			// microM_lkm as l=3
//  matrix microSigma;

};
 
} // namespace dem

#endif
