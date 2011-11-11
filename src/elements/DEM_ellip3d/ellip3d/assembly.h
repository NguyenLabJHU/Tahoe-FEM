#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "vec.h"
#include "gradation.h"
#include "particle.h"
#include "contact.h"
#include "boundary.h"
#include "rectangle.h"
#include "cylinder.h"
#include <map>
#include <list>
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
	    Volume = 0;
	    BulkDensity = 0;
	    PossCntctNum = 0;
	    ActualCntctNum = 0;
	    RgdBdryNum = 0;
	    FlbBdryNum = 0;
	    Gravity = false;
	};
	
	~assembly(){
		std::list<particle*>::iterator pt;
		std::list<RGDBDRY*>::iterator  rt;
		std::list<FLBBDRY*>::iterator  ft;

		// it is important to release memory pointed to by pointers in the container,
		// otherwise memory leaking occurs
		for(pt=ParticleList.begin();pt!=ParticleList.end();++pt)
			delete (*pt);
		for(rt=RBList.begin();rt!=RBList.end();++rt)
			delete (*rt);
		for(ft=FBList.begin();ft!=FBList.end();++ft)
			delete (*ft);
	};

	particle* getParticle(int id) const{
		std::list<particle*>::const_iterator it;
		// find() algorithm is a better choice
		for (it=ParticleList.begin();it!=ParticleList.end();++it)
			if ((*it)->getID()==id) return (*it);
		return NULL;
	}

	void        createSample(const char* str);           // create a sample with particles from an existing file
	void        createRigidBoundary(std::ifstream &ifs); // create rigid boundaries from an existing file
	void        createFlexiBoundary(std::ifstream &ifs); // create flxible boundaries from an existing file
	void        createBoundary(const char* str);         // create either rigid or flexible boundaries from an existing file
	void        findContact();                           // detect and resolve contact between particles
	void        findParticleOnBoundary();                // find particles on boundaries
	void        findParticleOnLine();                    // find particles on lines
	void        createFlbNet();

	void        clearForce();                            // clear forces and moments for all particles
	void        flexiBoundaryForceZero();
	void        initFBForce();
	void        internalForce(long double& avgnm, long double& avgsh); // calculate inter-particle forces
	void        rigidBoundaryForce();                    // calcualte forces between rigid boundaries and particles
	void        rigidBoundaryForce(long double penetr[],int cntnum[]);
	void        flexiBoundaryForce();
	void        updateParticle();                        // update motion of particles

	long double ellipPileForce();                        // for force pile only
	void        ellipPileUpdate();                       // for force pile only

	vec         ellipPileDimn();
	long double ellipPileTipZ();
	long double ellipPilePeneVol();


	// if bn[i]=2, the 2nd rigid boundary should be updated according to rbctl[i],
	// totally num rigid boundaries must be updated
	void        updateRB(int bn[], UPDATECTL rbctl[], int num);
	void        updateRB6();
	void        updateRectPile();
	
	// if bn[i]=2, the 2nd flxible boundary should be updated according to fbctl[bn[i]*2-2] and fbctl[bn[i]*2-1], 
	// totally num flxible boundaries must be updated, the size of fbctl is 2 times large as size of bn
	void        updateFB(int bn[], UPDATECTL fbctl[], int num);

	long double getDensity() const; 
	int         getPossCntctNum() const {return  PossCntctNum;};
	int         getActualCntctNum() const {return ActualCntctNum;}
	long double getAveragePenetration() const;
	long double getMinTimeStep() const;
	long double getAverageVelocity() const;
	long double getAverageForce() const;
	long double getAverageOmga() const;
	long double getAverageMoment() const;
	long double getParticleVolume() const;
	vec         getTopFreeParticlePosition() const;
	long double getTransEnergy() const;
	long double getRotatEnergy() const;
	long double getKinetEnergy() const;
	long double getPotenEnergy(long double ref) const;

	vec         getNormalForce(int bdry) const;       // get normal force acting on the bdry_th rigid boundary
	vec         getShearForce(int bdry) const;        // get shear force acting on the bdry_th rigid boundary
	long double getAvgNormal(int bdry) const;
	vec         getApt(int bdry) const;               // get a point on bdry_th rigid boundary
	vec         getDirc(int bdry) const;              // get the dirc of bdry_th rigid boundry
	long double getArea(int bdry) const;
	long double getAverageRigidPressure() const;
	void        setArea(int bdry,long double a);      // set the area of the bdry-th rigid boundary be a

	void        printParticle(const char* str) const; // print particles info into a disk file
	void        printContact(const char* str) const;  // print contacts information
	void        printBoundary(const char* str) const; // print rigid boundaries info to a disk file
	void        printRectPile(const char* str);       // append rectangular pile info into a disk file

	// create a specimen by depositing particles into rigid boundaries
	void deposit_RgdBdry(gradation& grad,
			     int   freetype,
			     int   total_steps,  
			     int   snapshots,
			     int   interval,
			     long double height,
			     const char* iniptclfile,   
			     const char* inibdryfile,
			     const char* particlefile, 
			     const char* contactfile,
			     const char* progressfile, 
			     const char* creparticle,
			     const char* creboundary,
			     const char* debugfile);

	// create a specimen by depositing particles into particle boundaries
	void deposit_PtclBdry(gradation& grad,
			      int   freetype,
			      long double rsize,
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
			    long double dimn         =0.05,              // dimension of particle-composed-boundary
			    long double rsize        =1.0,               // relative container size
			    const char* iniptclfile  ="dep_particle_end",// input file, initial particles
			    const char* particlefile ="scl_particle",    // output file, resulted particles, including snapshots 
			    const char* contactfile  ="scl_contact",     // output file, resulted contacts, including snapshots
			    const char* progressfile ="scl_progress",    // output file, statistical info
			    const char* debugfile    ="scl_debug");      // output file, debug info

	// generate particles in space for rigid boundaries
	void generate(gradation& grad, 
		      const char* str,
		      int freetype,
		      long double ht);

	// generate particles in space for particle boundaries
	void generate_p(gradation& grad,
			const char* str,
			int freetype,
			long double rsize,
			long double ht);

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

	// actual deposit function for particle boundaries
	void deposit_p(int         total_steps  =50000,             // total_steps
		       int         snapshots    =100,               // number of snapshots   
		       int         interval     =10,                // print interval 
		       long double dimn   =0.05,                    // dimension of particle-composed-boundary
		       long double rsize  =1.0,                     // relative container size
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

	void collapse(int   rors, 
		      int   total_steps,  
		      int   snapshots,
		      int   interval,
		      const char* iniptclfile,
		      const char* initboundary,
		      const char* particlefile,
		      const char* contactfile,
		      const char* progressfile,
		      const char* debugfile);
	
	void setBoundary(int rors,
		     int bdrynum,
		     long double dimn,
		     const char* boundaryfile);

	void trim(int rors,
		  const char* iniptclfile,
		  const char* inibdryfile,
		  const char* particlefile,
		  const char* boundaryfile);

	void TrimPtclBdryByHeight(double height,
				  const char* iniptclfile,
				  const char* particlefile);

        // Isotropically compress floating particles to a specific confining pressure, which is usually a low
        // value in order to create an intial status. Force boundaries are used. This process may be not 
        // physically true.
	void isotropic(int          total_steps  =100000,
		       int          snapshots    =100,
		       int          interval     =10,
		       long double  sigma        =1.0e+4,
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
		       long double  sigma_a       =1.0e+4,
		       long double  sigma_b       =1.0e+5,	
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
		       long double  sigma_values[],
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
		      long double  sigma_3        =1.0e+4,
		      long double  sigma_1        =1.0e+5,
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
		      long double  sigma_values[],
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
				 double       sigma        =5.0e+5,
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
		      long double  sigma_a      =1.0e+5,
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
		      long double  sigma_a      =3.0e+5,
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
			    long double dimn         =0.05,
			    long double rsize        =1.0,
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
			      long double dimn         =0.05,
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
				long double dimn         =0.05,
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
			     long double dimn         =0.05,
			     long double force        =1.0e+4,
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
			  long double  sigma_a       =1.0e+4,
			  long double  sigma_w       =1.0e+5,
			  long double  sigma_l       =1.0e+5,	
			  long double  sigma_h       =1.0e+5,	
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
	
	void soft_tric(long double _sigma3, long double _b,
		       const char* iniptclfile ="isotropic",
		       const char* boundaryfile="isobdry",
		       const char* responsefile="sftc",
		       const char* resultfile  ="sftcompressed",
		       const char* trackfile   ="tracksft");

	void earthPressure(long double pressure, bool IsPassive,
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

	void simpleShear(long double normal_pressure,long double _b,
                         const char* iniptclfile ="isotropic",
                         const char* boundaryfile="isobdry",
                         const char* responsefile="simpleshear",
                         const char* resultfile  ="simplesheared",
                         const char* trackfile   ="tracksimple");

	void dircShear(long double rate, long double roterate, long double stress,
		       const char* iniptclfile ="deposit",
		       const char* boundaryfile="rgdcube.data",
		       const char* responsefile="ssvds",
		       const char* resultfile  ="dircshear",
		       const char* trackfile   ="trackds");

private:
	// gravitation property
	bool Gravity;        

	// particles property
	int  TotalNum;                      // total number of particles
	int  PossCntctNum;                  // possible contact number based on spherical distances
	int  ActualCntctNum;                // actual contact number based on solution of 6th order equations
	std::list<particle*>  ParticleList; // a list of pointers, each pointing to a particle
	std::list<CONTACT>    ContactList;  // a list of contacts
	std::vector<cnttgt>   CntTgtVec;    // a vector to store tangential contact force and displacement
	gradation             Gradation;    // particles gradation
 
	// container property
	int  RORC;               // rectangular--1 or cylindrical--0
	cylinder S;              // S - cylinder specimen
	rectangle R;             // R - rectangle specimen
	long double Volume;      // volume of the specimen
	long double BulkDensity; // bulk density of specimen

	// boundary property
	int  BdryType;              // 0 - rigid boundaries; 1 - flxible boundaries
	int  RgdBdryNum;            // rigid boundary number
	int  FlbBdryNum;            // flxible boundary number
	std::list<RGDBDRY*> RBList; // a list of pointers, each pointing to a rigid boundary
	std::list<FLBBDRY*> FBList; // a list of pointers, each pointing to a flexible boundary
	std::map<int,std::vector<boundarytgt>  > BdryTgtMap; // a map to store particle-boundary contacts' tangential info
};

} // namespace dem ends

#endif
