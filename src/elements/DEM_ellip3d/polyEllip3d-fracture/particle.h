#ifndef PARTICLE_H
#define PARTICLE_H

#include "realtypes.h"
#include "vec.h"
#include "gradation.h"
#include "contact.h"
#include "boundary.h"
#include "rectangle.h"
#include "cylinder.h"
#include "boundarytgt.h"
#include "matrix.h"
#include <map>

namespace dem {

  class particle{
    friend class assembly;
    friend class contact<particle>;
    friend class flb_bdry<particle>;
  public:
    particle(int n, int type, vec center_geo, REAL r, REAL young, REAL poisson);	// August 16, 2013
    particle(int n, int type, vec center_geo, REAL aplus, REAL aminus, REAL bplus, REAL bminus, REAL cplus, REAL cminus, REAL young, REAL poisson);	// August 16, 2013
    particle(int id,int type, REAL aplus, REAL aminus, REAL bplus, REAL bminus, REAL cplus, REAL cminus, vec position_geo, vec dirca, vec dircb, vec dircc, REAL young, REAL poisson);	// August 16, 2013
    particle(int n, int type, vec center, gradation& grad, REAL young, REAL poisson);	// August 19, 2013
    particle(const particle& pt, int brk_pln);	// September 23, 2014. for fracture model 
    
    int  getID() const {return ID;}
    int  getType() const {return type;}
    int  getCandidacy() const {return candidacy;}	// for particle fracture, October 17, 2013
    int  getNumBroken() const {return numBroken;}
    int  getNumCriticalContacts() const {return numCriticalContacts;}
    int  getTypeBroken() const {return typeBroken;}
    int  getBreakPlane();	// determine which plane is break plane by average stress, September 22, 2014 
    int  getBreakPlane(int abc);	// determine which plane is break plane by average stress, September 22, 2014 
					// int abc is one of the two axles that construct the plane, for case (5)
    int  calculateBreakPlane();	// this will first calculate maximum tensile stress
				// then return break plane [-1, 1,2,3], -1 means not break this particle
    void rotatePrincipalDirections(REAL V[3][3]);// for spheres, rotate their geometry principal directions the same as the
					// principal directions of the stress
    void rotatePrincipalDirections(dem::vec);// for spheres, rotate their geometry principal directions the same as the
					// break plane determined by the three maximum tensile stress in contacts
    void rotateToXYZDirections();	// rotate the particles to be along x,y,z directions
    void rotateBreakPlaneRandomly();	// rotate break plane randomly a little
    dem::vec calculateInitialCohesiveForce();	// calculate the initial cohesive force provided by 4 springs on particle p1
						// based on the fact that when particle is sub-divided, then two sub-poly-ellipsoids
						// should keep the same acceleration
    REAL getAplus() const {return aplus;}
    REAL getAminus() const {return aminus;}
    REAL getBplus() const {return bplus;}	// August 16, 2013
    REAL getBminus() const {return bminus;}
    REAL getCplus() const {return cplus;}
    REAL getCminus() const {return cminus;}
    REAL getAverageRadius() const {return (aplus+aminus+bplus+bminus+cplus+cminus)*0.166666666666667;}
    REAL getAMax() const {return aplus>aminus?aplus:aminus;}
    REAL getBMax() const {return bplus>bminus?bplus:bminus;}
    REAL getCMax() const {return cplus>cminus?cplus:cminus;}
    REAL getAMin() const {return aplus>aminus?aminus:aplus;}
    REAL getBMin() const {return bplus>bminus?bminus:bplus;}
    REAL getCMin() const {return cplus>cminus?cminus:cplus;}
    REAL getMaxRadius() const;	// get max among aplus, aminus, bplus, bminus, cplus and cminus. August 21, 2013
    REAL getMinRadius() const;	// get min. August 27, 2013
    REAL getYoung() const {return young;}
    REAL getPoisson() const {return poisson;};
    REAL getVolume() const {return volume;}
    REAL getMass() const {return mass;}
    REAL getDensity() const {return density;}
    REAL getStress1() const {return stress1;}
    REAL getStress2() const {return stress2;}
    REAL getStress3() const {return stress3;}
    REAL getStrengthContact() const {return strengthContact;}
    vec  getCurrPosition() const {return curr_position;}
    vec  getPrevPosition() const {return prev_position;}
    vec  getInitCenterMass() const {return init_center_mass;}	// initial position for granular strain
    vec  getStartCenterMass() const {return start_center_mass;}
    vec  getCurrCenterMass() const {return curr_center_mass;}	// August 19, 2013
    vec  getPrevCenterMass() const {return prev_center_mass;}
    vec  getCurrDirecA() const {return curr_direction_a;}
    vec  getCurrDirecB() const {return curr_direction_b;}
    vec  getCurrDirecC() const {return curr_direction_c;}
    vec  getPrevDirecA() const {return prev_direction_a;}
    vec  getPrevDirecB() const {return prev_direction_b;}
    vec  getPrevDirecC() const {return prev_direction_c;}
    vec  getCurrVelocity() const {return curr_velocity;}
    vec  getPrevVelocity() const {return prev_velocity;}
    vec  getCurrOmga() const {return curr_omga;}
    vec  getPrevOmga() const {return prev_omga;}
    vec  getCurrAcceleration() const {return curr_acceleration;}
    vec  getPrevAcceleration() const {return prev_acceleration;}
    vec  getForce() const {return force;}
    vec  getMoment() const {return moment;}
    vec  getConstForce() const {return const_force;}
    vec  getConstMoment() const {return const_moment;}
    vec  getJ() const {return J;}
    bool isInContact() const {return inContact;}
    bool getIsAbleDivide() const {return isAbleDivide;}
    bool isCLongEnough() const;	// if c^+ + c^- > 1.2max(A_max, B_max)
    bool isBLongEnough() const;	// if b^+ + b^- > 1.2max(A_max, C_max)
    bool isALongEnough() const;	// if a^+ + a^- > 1.2max(B_max, C_max)
    bool isBLongerThanC() const {return bplus+bminus>1.2*(cplus+cminus);}
    bool isCLongerThanB() const {return cplus+cminus>1.2*(bplus+bminus);}
    bool isALongerThanC() const {return aplus+aminus>1.2*(cplus+cminus);}
    bool isCLongerThanA() const {return cplus+cminus>1.2*(aplus+aminus);}
    bool isALongerThanB() const {return aplus+aminus>1.2*(bplus+bminus);}
    bool isBLongerThanA() const {return bplus+bminus>1.2*(aplus+aminus);}

    REAL getRadius(vec v, int, int, int) const;	// August 21, 2013
//    REAL getRadius(vec v);
    REAL getTransEnergy() const;
    REAL getRotatEnergy() const;
    REAL getKinetEnergy() const;
    REAL getPotenEnergy(REAL ref) const;	// August 16, 2013
    
    void setID(int n){ID=n;}
    void setType(int n) {type=n;}
    void candidacyMinus() {candidacy--;}	// for particle fracture, October 17, 2013
    void candidacyPlus() {candidacy++;} 
    void setAplus(REAL dd){aplus=dd;}	
    void setAminus(REAL dd){aminus=dd;}
    void setBplus(REAL dd){bplus=dd;}	// August 16, 2013
    void setBminus(REAL dd){bminus=dd;}
    void setCplus(REAL dd){cplus=dd;}
    void setCminus(REAL dd){cminus=dd;}
    void setIsAbleDivideFalse(){isAbleDivide = false;}
    void expand(REAL percent) {aplus *= (1+percent); aminus *= (1+percent); bplus *= (1+percent); bminus *= (1+percent); cplus *= (1+percent); cminus *= (1+percent);}	// August 16, 2013
    void setCurrPosition(vec vv){curr_position=vv;}
    void setPrevPosition(vec vv){prev_position=vv;}
    void setCurrCenterMass(vec vv){curr_center_mass=vv;}	// August 19, 2013
    void setPrevCenterMass(vec vv){prev_center_mass=vv;}
    void setInitCenterMass(vec vv){init_center_mass=vv;}	// initial center for granular strain
    void setStartCenterMass(vec vv){start_center_mass=vv;}	// position at starting step when tessellation is regenerated
    void setCurrDirecA(vec vv){curr_direction_a=vv;}
    void setCurrDirecB(vec vv){curr_direction_b=vv;}
    void setCurrDirecC(vec vv){curr_direction_c=vv;}
    void setPrevDirecA(vec vv){prev_direction_a=vv;}
    void setPrevDirecB(vec vv){prev_direction_b=vv;}
    void setPrevDirecC(vec vv){prev_direction_c=vv;}
    void setCurrVelocity(vec vv){curr_velocity=vv;}
    void setPrevVelocity(vec vv){prev_velocity=vv;}
    void setCurrOmga(vec vv){curr_omga=vv;}
    void setPrevOmga(vec vv){prev_omga=vv;}
    void setCurrAcceleration(vec vv){curr_acceleration=vv;}
    void setPrevAcceleration(vec vv){prev_acceleration=vv;}
    void setForce(vec vv){force=vv;}
    void setMoment(vec vv){moment=vv;}
    void setConstForce(vec vv){const_force=vv;}
    void setConstMoment(vec vv){const_moment=vv;}
    void setJ(vec v){J=v;}
    void setMass(REAL d){mass=d;}
    void setDensity(REAL dn) {density=dn;}
    void setExternForce(vec fc) {force=fc;}
    void setExternMoment(vec mm) {moment=mm;} 
    void initTypeBroken() {typeBroken=0;}
    
    void clearForce();	// Add moment by gravity. August 28, 2013
    void addForce(vec vv) {force+=vv;}	// April 23, 2014
    void addMoment(vec vv) {moment+=vv;}
    void clearStress();
    void addStress(matrix mm) {average_stress = average_stress+mm;}
    void calcStress();	// calculate the average stress after the summation of force terms
    matrix getStress() {return average_stress;}
    void addMaximuContactTensile(REAL, dem::vec);
    void update();	// August 16, 2013
    void shapeTransition();	// transit bad-shaped poly-ellipsoids generated in fracture
				// to normal shaped poly-ellipsoids, suggested by Dr. Fu at LLNL  
    void calculateGeometry();	// calculate the geometry property of the new transited poly-ellipsoid 

    // update global coefficients in the following form based on position/dimensions/orientations
    // a0 x^2 + a1 y^2 + a2 z^2 + a3 xy + a4 yz + a5 zx + a6 x + a7 y + a8 z + a9 = 0
//    void GlobCoef(int, int, int);  // August 21, 2013
    void GlobCoef();	// September 6, 2013. Calculate coefficients for all 8 octants at every time step
    void getGlobCoef(REAL coef[], int num_oct) const; // fetch global coeffs into coef[].. September 6, 2013
    REAL surfaceError(vec pt) const;
    vec  localVec(vec) const;   // transform a vector in global coordinates into local coordinates
    vec  globalVec(vec) const;  // transform a vector in local coordinates into global coordinates
    void print() const;	// August 16, 2013
    
    // Assumption: a particle only intersects a plane a little and it cannot pass through the plane
    //             with most of its body, this is guaranteed by contacting forces.
    
    //v is the point the line passing through
    //d is the unit vector parallel to the line
    bool intersectWithLine(vec v, vec dirc, vec rt[], int num_oct) const;
    
    //find the point on plane which is deepest into a particles, px+qy+rz+s=0 is the equation of the plane
    //true means intersection; false means no intersection.
    bool nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, vec& ptnp, int num_oct) const;	// August 19, 2013
    
    //side indicates which side the particles are in about the plane
    //calculate the normal force between particle and a plane rigid boundary
    void planeRBForce(plnrgd_bdry<particle>* plb,	// August 19, 2013. apply moment on the mass center. August 22, 2013
		      std::map<int,std::vector<boundarytgt> >& BoundarytgtMap,
		      std::vector<boundarytgt>& vtmp,
		      REAL &penetr);
    
    //if side<0, particles are inside the cylinder; else side is outside the cylinder
    //calculate the normal force between particle and a cylinder wall
    vec cylinderRBForce(int bdry_id, const cylinder& S, int side);	// August 19, 2013. August 22, 2013

    void breakItSelf(int break_plane);	// change the particle itself to the fractured particle which is breaked by the break_plane of the original particle
					// and also in the positive part, September 22, 2014. 
    
  private:
    // particle type (settings for individual particle):
    //  0 - free particle
    //  1 - fixed particle
    //  2 - special case 2 (pure moment): translate first, then rotate only, MNT_START needs to be defined
    //  3 - special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
    //  4 - special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only
    //  5 - free boundary particle
    // 10 - ghost particle
    int ID;
    int type;            
    REAL aplus, aminus, bplus, bminus, cplus, cminus;   // six parameters, not order by value, corresponding to three principal axels, August 16, 2013.
    REAL young;
    REAL poisson;
    vec curr_position;   // particle center, center of geometry August 16, 2013
    vec prev_position;
    vec local_center_mass;	// local coordinate of mass center with center_geo as origin and princinpal directions as axels, unchanged. August 16, 2013
    vec prev_center_mass;
    vec curr_center_mass;	//global coordinates of previous and current center of mass
    vec init_center_mass;	// initial center
    vec start_center_mass;	// position at starting step when tessellation is regenerated
    vec curr_direction_a, curr_direction_b, curr_direction_c;// the direction of the three axles, in radian. it is the angle, not directional vector
    vec prev_direction_a, prev_direction_b, prev_direction_c;
    vec curr_velocity;   // the velocity of the mass center
    vec prev_velocity;
    vec curr_omga;       // angular velocity in global frame!
    vec prev_omga;
    vec curr_acceleration;
    vec prev_acceleration;
    vec curr_acce_rotate;	// rotation acceleration
    vec prev_acce_rotate;
    vec pre_force;
    vec force;
    vec pre_moment;
    vec moment;
    vec const_force;
    vec const_moment;
    vec flb_force;
    vec flb_moment;
    vec mres;            // resistence moment provided by normal distribution force
    matrix average_stress;	// the particle average stress based on the dynamic version of Andrade's formula
    REAL density; // specific gravity
    REAL mass;
    REAL volume;
    vec  J;              // moment of inertia in local body-fixed frame
    REAL coef1[10];// record particle's coefficients in global coordinates
    REAL coef2[10];
    REAL coef3[10];
    REAL coef4[10];
    REAL coef5[10];	// September 6, 2013
    REAL coef6[10];
    REAL coef7[10];
    REAL coef8[10];
    REAL kinetEnergy; // kinetic energy
    int  cntnum;
    bool inContact; // in contact with other particle or boundary
    int candidacy; // for particle fracture, October 17, 2013. Implies the number of fracture plane of this particle. Can be more than 1 if allow further-level subdivision.
		   // now this value is only 1 or 0 in the new fracture model 
    void init();	// August 19, 2013   
    bool isAbleDivide;	// the particle can be divided or not, particles of case(2) and case(4) cannot be divided, September 22, 2014

    // for the fracture criterion based on contact tensile stress
    int numCriticalContacts;	// the number of critical contact points, which is the contact whose maximum tensile stress
				// is larger than the critical tensile strength.
    REAL stress1;	// the maximum tensile stress for the three contact points
    REAL stress2;	// whose maximum tensile stresses are larger than the tensile strength
    REAL stress3;	// if 0, means this is not a critical contact point
    dem::vec contact1;	// the global coordinate of the contact point corresponding to stress1
    dem::vec contact2;	// has and only has three contact points, it is because a plane can be determined by 3 points
    dem::vec contact3;  // Even we have 4th contact points who is also critical point, we will not account for it

    int numBroken;	// the number of brokens that this particle has experienced
    int typeBroken;	// type of this broken, -1 means Hoek-Brown broken, 0 means not break, 1 means 3 contact points broken
			// 2 means 2 contact points broken

    // the movement increment for the bad-shaped poly-ellipsoid in the transition process
    // currently, we only transit these fracture particles
    REAL delta_a;	// the increment in b axle, initially it is 0, means not transit
    REAL delta_b;	// positive means, move in the positive direction in local system
    REAL delta_c;	// negative means, move in the negative direction in local system

    REAL weibullPhi;	// random number between 0 and 1, defined when the particle is generated
			// represent the weibull statistics of strength
    REAL strengthHoek;		// strength for Hoek-Brown criterion
    REAL strengthContact;	// strength for the maximum tensile stress in contacts criterion

    // the integrals over poly-ellipsoids for average stress calculation
    REAL int_x, int_y, int_z, int_xy, int_xz, int_yz, int_x2, int_y2, int_z2;

  public:
    int IsFBP;
  };
  
} // namespace dem ends

#endif
