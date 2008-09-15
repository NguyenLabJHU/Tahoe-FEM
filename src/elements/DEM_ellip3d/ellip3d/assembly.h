#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <map>
#include <list>
#include <vector>
#include <fstream>

#include "vec.h"
#include "gradation.h"
#include "particle.h"
#include "contact.h"
#include "boundary.h"
#include "rectangle.h"
#include "cylinder.h"

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

		// it is essential to release memory pointed to by pointers in the container,
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


	void createSample(char* str);
	void createContact();
	void createRB(std::ifstream &ifs); // create rigid boundary from an existing file
	void createFB(std::ifstream &ifs); // create flxible boundary from an existing file
	void createBdry(char* str);        // create both flxible boundary and rigid boundary from an existing file
	void createPBL();
	void createPLL();
	void createFlbNet();

	void setForceZero();
	void fbForceZero();
	void initFBForce();
	void internForce(long double& avgnm, long double& avgsh);
	void rbForce();
	void rbForce(long double penetr[],int cntnum[]);
	void fbForce();
	void particleUpdate();

	long double ellipPileForce();  // for force pile only
	void        ellipPileUpdate(); // for force pile only

	vec         ellipPileDimn();
	long double ellipPileTipZ();
	long double ellipPilePeneVol();


	// if bn[i]=2, the 2nd rigid boundary should be updated according to rbctl[i],
	// totally num rigid boundaries must be updated
	void updateRB(int bn[], UPDATECTL rbctl[], int num);
	void updateRB6();
	void updateRectPile();
	
	// if bn[i]=2, the 2nd flxible boundary should be updated according to fbctl[bn[i]*2-2] and fbctl[bn[i]*2-1], 
	// totally num flxible boundaries must be updated, the size of fbctl is 2 times large as size of bn
	void updateFB(int bn[], UPDATECTL fbctl[], int num);

	long double density() const; //get density of the specimen in kg/m^3
	int  getPossCntctNum() const {return  PossCntctNum;};
	int  getActualCntctNum() const {return ActualCntctNum;}
	long double avgPenetration() const;
	long double avgVelocity() const;
	long double avgForce() const;
	long double avgOmga() const;
	long double avgMoment() const;
	long double ptclVolume() const;
	vec         topFreePtclPos() const;
	long double transEnergy() const;
	long double rotaEnergy() const;
	long double kinEnergy() const;
	long double potEnergy(long double ref) const;

	vec  getNormal(int bdry) const;          // get the pressure acting on the bdry_th rigid boundary
	vec  getTangt(int bdry) const;           // get the shear friction of the rigid boundary
	long double getAvgNormal(int bdry) const;
	vec  getApt(int bdry) const;             // get a point on bdry_th rigid boundary
	vec  getDirc(int bdry) const;            // get the dirc of bdry_th rigid boundry
	long double getArea(int bdry) const;
	long double avgRgdPres() const;
	void setArea(int bdry,long double a);    // set the area of the bdry-th rigid boundary be a

	void snapshot(char* str) const;          // print particles dynamic info into a disk file
	void printPtcl(char* str) const;         // print particles info into a disk file
	void printCntct(char* str) const;        // print contacts information
	void printBdry(char* str) const;         // print rigid boundaries info to a disk file
	void printRectPile(char* str);           // append rectangular pile info into a disk file
	void dispBdry() const;                   // display both rigid and flexible boundaries information

	// create a specimen from discreate particles through floating and then gravitation,
	// file cre_particle contains the final particle information,
	// file cre_boundary contains the final boundary information.
	void deposit_RgdBdry(gradation& grad,
			     int   freetype,
			     int   total_steps,  
			     int   snapshots,
			     long double height,
			     char* iniptclfile,   
			     char* inibdryfile,
			     char* particlefile, 
			     char* contactfile,
			     char* progressfile, 
			     char* creparticle,
			     char* creboundary,
			     char* exceptionfile);

	// create a specimen from discreate particles through floating and then gravitation,
	// but the boundaries are composed of fixed particles.
	// file cre_particle contains the final particle information,
	// file cre_boundary contains the final boundary information.
	void deposit_PtclBdry(gradation& grad,
			      int   freetype,
			      long double rsize,
			      int   total_steps,  
			      int   snapshots,
			      char* iniptclfile,   
			      char* particlefile, 
			      char* contactfile,
			      char* progressfile, 
			      char* exceptionfile);
	
	// scaling the whole assembly from deposited state, until it reaches steady state
	void scale_PtclBdry(int   total_steps  =50000,             // total_steps
			    int   snapshots    =100,               // number of snapshots   
			    long double dimn   =0.05,              // dimension of particle-composed-boundary
			    long double rsize  =1.0,               // relative container size
			    char* iniptclfile  ="dep_particle_end",// input file, initial particles
			    char* particlefile ="scl_particle",    // output file, resulted particles, including snapshots 
			    char* contactfile  ="scl_contact",     // output file, resulted contacts, including snapshots
			    char* progressfile ="scl_progress",    // output file, progress statistic information
			    char* exceptionfile="scl_exception");  // output file, progress float exceptions

	void init(gradation& grad, 
		  char* str,
		  int freetype,
		  long double ht);

	void init_p(gradation& grad,
		    char* str,
		    int freetype,
		    long double rsize,
		    long double ht);

        // actual deposit function for the case of rigid boundaries
	// the container can be as simple as a bottom plate
	void deposit(int   total_steps  =100000,            // total_steps
		     int   snapshots    =100,               // number of snapshots   
		     char* iniptclfile  ="flo_particle_end",// input file, initial particles
		     char* inibdryfile  ="dep_boundary_ini",// input file, initial boundaries
		     char* particlefile ="dep_particle",    // output file, resulted particles, including snapshots 
		     char* contactfile  ="dep_contact",     // output file, resulted contacts, including snapshots
		     char* progressfile ="dep_progress",    // output file, progress statistic information
		     char* exceptionfile="dep_exception");  // output file, progress float exceptions

	// actual deposit function for the case of fixed particle boundaries
	void deposit_p(int   total_steps  =50000,             // total_steps
		       int   snapshots    =100,               // number of snapshots   
		       long double dimn   =0.05,              // dimension of particle-composed-boundary
		       long double rsize  =1.0,               // relative container size
		       char* iniptclfile  ="flo_particle_end",// input file, initial particles
		       char* particlefile ="dep_particle",    // output file, resulted particles, including snapshots 
		       char* contactfile  ="dep_contact",     // output file, resulted contacts, including snapshots
		       char* progressfile ="dep_progress",    // output file, progress statistic information
		       char* exceptionfile="dep_exception");  // output file, progress float exceptions

        //squeeze paticles inside a container by moving the boundaries
	void squeeze(int   total_steps  =20000,             // total_steps
		     int   init_steps   =5000,              // initial_steps to reach equilibrium
		     int   snapshots    =100,               // number of snapshots   
		     int   flag         =-1,                // -1 squeeze; +1 loosen
		     char* iniptclfile  ="flo_particle_end",// input file, initial particles
		     char* inibdryfile  ="dep_boundary_ini",// input file, initial boundaries
		     char* particlefile ="dep_particle",    // output file, resulted particles, including snapshots 
		     char* boundaryfile ="dep_boundary",    // output file, resulted boundaries
		     char* contactfile  ="dep_contact",     // output file, resulted contacts, including snapshots
		     char* progressfile ="dep_progress",    // output file, progress statistic information
		     char* exceptionfile="dep_exception");  // output file, progress float exceptions

	void collapse(int   rors, 
		      int   total_steps,  
		      int   snapshots,
		      char* iniptclfile,
		      char* initboundary,
		      char* particlefile,
		      char* contactfile,
		      char* progressfile,
		      char* exceptionfile);
	
	void setbdry(int rors,
		     int bdrynum,
		     long double dimn,
		     char* boundaryfile);

	void trim(int rors,
		  char* iniptclfile,
		  char* inibdryfile,
		  char* particlefile,
		  char* boundaryfile);

	void TrimPtclBdryByHeight(double height,
				  char* iniptclfile,
				  char* particlefile);

        // Isotropically compress floating particles to a specific confining pressure, which is usually a low
        // value in order to create an intial status. Force boundaries are used. This process may be not 
        // physically true.
	void isotropic(int    total_steps  =100000,
		       int    snapshots    =100,
		       long double sigma   =1.0e+4,
		       char*  iniptclfile  ="flo_particle_end",
		       char*  inibdryfile  ="iso_inbdry",
		       char*  particlefile ="iso_particle",
		       char*  boundaryfile ="iso_boundary",
		       char*  contactfile  ="iso_contact",
		       char*  progressfile ="iso_progress",
		       char*  balancedfile ="iso_balanced",
		       char*  exceptionfile="iso_exception");

        // The specimen has been isotropically compressed to confining pressure sigma_a. This function
        // increases confining pressure step by step to sigma_b, thus making it possible to find out
        // balanced status where particle pressure equals confining pressure. Force boundaries are used.
	void isotropic(int    total_steps   =100000,
		       int    snapshots     =100,
		       long double sigma_a  =1.0e+4,
		       long double sigma_b  =1.0e+5,	
		       int    sigma_division=100,	  
		       char*  iniptclfile   ="iso_particle_10k",
		       char*  inibdryfile   ="iso_boundary_10k",
		       char*  particlefile  ="iso_particle", 
		       char*  boundaryfile  ="iso_boundary", 
		       char*  contactfile   ="iso_contact",
		       char*  progressfile  ="iso_progress",
		       char*  balancedfile  ="iso_balanced", 
		       char*  exceptionfile ="iso_exception");
	
        // The specimen has been isotropically compressed to confining pressure sigma_a. This function
	// follows an unloading-reloading stress path. Force boundaries are used.
	void isotropic(int    total_steps,
		       int    snapshots,
		       int    sigma_points,			  
		       long double sigma_values[],
		       int    sigma_division=100,
		       char*  iniptclfile   ="iso_particle_10k",
		       char*  inibdryfile   ="iso_boundary_10k",
		       char*  particlefile  ="iso_particle", 
		       char*  boundaryfile  ="iso_boundary", 
		       char*  contactfile   ="iso_contact",
		       char*  progressfile  ="iso_progress",
		       char*  balancedfile  ="iso_balanced", 
		       char*  exceptionfile ="iso_exception");
	
        // The specimen has been isotropically compressed to confining pressure sigma_3. This function
        // increases confining pressure step by step to sigma_1, thus making it possible to find out
        // balanced status where top & bottom particle pressure equals major principle stress. 
        // Side boundaries are fixed, top and bottom plates are force-controlled.
	void odometer(int    total_steps    =100000,
		      int    snapshots      =100,
		      long double sigma_3   =1.0e+4,
		      long double sigma_1   =1.0e+5,
		      int    sigma_division =100,		  
		      char*  iniptclfile    ="iso_particle_10k",
		      char*  inibdryfile    ="iso_boundary_10k",
		      char*  particlefile   ="odo_particle", 
		      char*  boundaryfile   ="odo_boundary", 
		      char*  contactfile    ="odo_contact",
		      char*  progressfile   ="odo_progress",
		      char*  balancedfile   ="odo_balanced", 
		      char*  exceptionfile  ="odo_exception");

        // The specimen has been isotropically compressed to confining pressure sigma_3. This function
        // increases confining pressure step by step to sigma_1, thus making it possible to find out
        // balanced status where top & bottom particle pressure equals major principle stress. 
        // Side boundaries are fixed, top and bottom plates are force-controlled. Unloading is applied.
	void odometer(int    total_steps,
		      int    snapshots,
		      int    sigma_points,			  
		      long double sigma_values[],
		      int    sigma_division=100,		  
		      char*  iniptclfile   ="iso_particle_10k",
		      char*  inibdryfile   ="iso_boundary_10k",
		      char*  particlefile  ="odo_particle", 
		      char*  boundaryfile  ="odo_boundary", 
		      char*  contactfile   ="odo_contact",
		      char*  progressfile  ="odo_progress",
		      char*  balancedfile  ="odo_balanced", 
		      char*  exceptionfile ="odo_exception");

        // The confining pressure is 500kPa. This function initializes triaxial compression test.
	void triaxialPtclBdryIni(int    total_steps  =10000,
				 int    snapshots    =100,
				 double sigma        =5.0e+5,
				 char*  iniptclfile  ="ini_particle_ini",
				 char*  inibdryfile  ="ini_boundary_ini",
				 char*  particlefile ="ini_particle", 
				 char*  boundaryfile ="ini_boundary", 
				 char*  contactfile  ="ini_contact",
				 char*  progressfile ="ini_progress",
				 char*  exceptionfile="ini_exception");

        // The confining pressure is 500kPa. This function performs triaxial compression test.
        // Displacement boundaries are used in axial direction.
	void triaxialPtclBdry(int    total_steps  =100000,
			      int    snapshots    =100,
			      char*  iniptclfile  ="iso_particle_100k",
			      char*  inibdryfile  ="iso_boundary_100k",
			      char*  particlefile ="tri_particle", 
			      char*  boundaryfile ="tri_boundary", 
			      char*  contactfile  ="tri_contact",
			      char*  progressfile ="tri_progress",
			      char*  balancedfile ="tri_balanced", 
			      char*  exceptionfile="tri_exception");

        // The specimen has been isotropically compressed to confining pressure sigma_a. This function
        // performs triaxial compression test. Displacement boundaries are used in axial direction.
	void triaxial(int    total_steps  =100000,
		      int    snapshots    =100,
		      long double sigma_a =1.0e+5,
		      char*  iniptclfile  ="iso_particle_100k",
		      char*  inibdryfile  ="iso_boundary_100k",
		      char*  particlefile ="tri_particle", 
		      char*  boundaryfile ="tri_boundary", 
		      char*  contactfile  ="tri_contact",
		      char*  progressfile ="tri_progress",
		      char*  balancedfile ="tri_balanced", 
		      char*  exceptionfile="tri_exception");
	
        // The specimen has been isotropically compressed to confining pressure sigma_a. This function
        // performs triaxial compression test with unloading. Displacement boundaries are used in 
        // axial direction.
	void triaxial(int    total_steps  =200000,
		      int    unload_step  =100000,
		      int    snapshots    =100,
		      long double sigma_a =3.0e+5,
		      char*  iniptclfile  ="iso_particle_300k",
		      char*  inibdryfile  ="iso_boundary_300k",
		      char*  particlefile ="tri_particle", 
		      char*  boundaryfile ="tri_boundary", 
		      char*  contactfile  ="tri_contact",
		      char*  progressfile ="tri_progress",
		      char*  balancedfile ="tri_balanced", 
		      char*  exceptionfile="tri_exception");
	
        // The specimen has been deposited with gravitation within boundaries composed of particles.
        // A rectangular pile is then drived into the particles using displacement control.
	void rectPile_Disp(int    total_steps  =50000,
			   int    snapshots    =100,
			   char*  iniptclfile  ="pile_particle_ini",
			   char*  inibdryfile  ="pile_boundary_ini",
			   char*  particlefile ="pile_particle", 
			   char*  boundaryfile ="pile_boundary", 
			   char*  contactfile  ="pile_contact",
			   char*  progressfile ="pile_progress",
			   char*  exceptionfile="pile_exception");
	
        // The specimen has been deposited with gravitation within boundaries composed of particles.
        // An ellipsoidal pile is then drived into the particles using displacement control.
	void ellipPile_Disp(int   total_steps  =50000,  
			    int   snapshots    =100, 
			    long double dimn   =0.05,
			    long double rsize  =1.0,
			    char* iniptclfile  ="pile_particle_ini",
			    char* particlefile ="pile_particle", 
			    char* contactfile  ="pile_contact",  
			    char* progressfile ="pile_progress",
			    char* exceptionfile="pile_exception" );

        // The specimen has been deposited with gravitation within rigid boundaries.
        // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
	void ellipPile_Impact(int   total_steps  =50000,  
			      int   snapshots    =100, 
			      long double dimn   =0.05,
			      char* iniptclfile  ="ipt_particle_ini",
			      char* inibdryfile  ="dep_boundary_ini",
			      char* particlefile ="ipt_particle", 
			      char* contactfile  ="ipt_contact",  
			      char* progressfile ="ipt_progress",
			      char* exceptionfile="ipt_exception" );

        // The specimen has been deposited with gravitation within particle boundaries.
        // An ellipsoidal penetrator is then impacted into the particles with initial velocity.
	void ellipPile_Impact_p(int   total_steps  =50000,  
				int   snapshots    =100, 
				long double dimn   =0.05,
				char* iniptclfile  ="ipt_particle_ini",
				char* particlefile ="ipt_particle", 
				char* contactfile  ="ipt_contact",  
				char* progressfile ="ipt_progress",
				char* exceptionfile="ipt_exception" );

        // The specimen has been deposited with gravitation within boundaries composed of particles.
        // An ellipsoidal pile is then drived into the particles using force control.
	void ellipPile_Force(int   total_steps  =50000,  
			     int   snapshots    =100, 
			     long double dimn   =0.05,
			     long double force  =1.0e+4,
			     int   division     =100,
			     char* iniptclfile  ="pile_particle_ini",
			     char* particlefile ="pile_particle", 
			     char* contactfile  ="pile_contact",  
			     char* progressfile ="pile_progress",
			     char* balancedfile ="pile_balanced",
			     char* exceptionfile="pile_exception" );

	void truetriaxial(int    total_steps   =1000000,
			  int    snapshots     =100,
			  long double sigma_a  =1.0e+4,
			  long double sigma_w  =1.0e+5,
			  long double sigma_l  =1.0e+5,	
			  long double sigma_h  =1.0e+5,	
			  int    sigma_division=100,			  
			  char*  iniptclfile   ="iso_particle_10k",
			  char*  inibdryfile   ="iso_boundary_10k",
			  char*  particlefile  ="tru_particle", 
			  char*  boundaryfile  ="tru_boundary", 
			  char*  contactfile   ="tru_contact",
			  char*  progressfile  ="tru_progress",
			  char*  balancedfile  ="tru_balanced", 
			  char*  exceptionfile ="tru_exception");

	void unconfined(int    total_steps  =100000,
			int    snapshots    =100,
			char*  iniptclfile  ="flo_particle_end",
			char*  inibdryfile  ="unc_inbdry",
		 	char*  particlefile ="unc_particle", 
			char*  contactfile  ="unc_contact",
			char*  progressfile ="unc_progress",
			char*  exceptionfile="unc_exception");
	
	void soft_tric(long double _sigma3, long double _b,
		       char* iniptclfile ="isotropic",
		       char* boundaryfile="isobdry",
		       char* responsefile="sftc",
		       char* resultfile  ="sftcompressed",
		       char* trackfile   ="tracksft");

	void earthPressure(long double pressure, bool IsPassive,
			   char* iniptclfile ="isotropic",
			   char* boundaryfile="isobdry",
			   char* responsefile="ssvpp",
			   char* resultfile  ="pp",
			   char* trackfile   ="pptrack");

	void shallowFoundation(char* iniptclfile ="isotropic",
			       char* boundaryfile="isobdry",
			       char* responsefile="ssvsf",
			       char* resultfile  ="shallowcompressed",
			       char* trackfile   ="shallowtrack");

	void simpleShear(long double normal_pressure,long double _b,
                         char* iniptclfile ="isotropic",
                         char* boundaryfile="isobdry",
                         char* responsefile="simpleshear",
                         char* resultfile  ="simplesheared",
                         char* trackfile   ="tracksimple");

	void dircShear(long double rate, long double roterate, long double stress,
		       char* iniptclfile ="deposit",
		       char* boundaryfile="rgdcube.data",
		       char* responsefile="ssvds",
		       char* resultfile  ="dircshear",
		       char* trackfile   ="trackds");

private:
	// gravitation property
	bool Gravity;        

	// particles property
	int  TotalNum;       // total number of particles
	int  PossCntctNum;   // possible contact number based on spherical distances
	int  ActualCntctNum; // actual contact number based on solution of 6th order equations
	std::list<particle*>  ParticleList;     // a list of pointers, each pointing to a particle
	std::list<CONTACT>    ContactList;      // a list of contacts
	std::vector<cnttgt>   CntTgtVec;      // a vector to store contacts' tangential force and displacement
	gradation             Gradation;        // particles gradation
 
	// container property
	int  RORC;              // rectangular--1 or cylindrical--0?
	cylinder S;             // S - cylinder specimen
	rectangle R;            // R - rectangle specimen
	long double Volume;     // volume of the specimen
	long double BulkDensity;// bulk density of specimen

	// boundary property
	int  BdryType;          // 0 - rigid boundaries; 1 - flxible boundaries.
	int  RgdBdryNum;        // rigid boundary number
	int  FlbBdryNum;        // flxible boundary number
	std::list<RGDBDRY*> RBList; // a list of pointers, each pointing to a rigid boundary
	std::list<FLBBDRY*> FBList; // a list of pointers, each pointing to a flexible boundary
	std::map<int,std::vector<boundarytgt>  > BdryTgtMap; // a map to store particle-boundary contacts' tangential info
};

} // namespace dem ends

#endif
