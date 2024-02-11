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
    
    int  getID() const {return ID;}
    int  getType() const {return type;}
    REAL getAplus() const {return aplus;}
    REAL getAminus() const {return aminus;}
    REAL getBplus() const {return bplus;}	// August 16, 2013
    REAL getBminus() const {return bminus;}
    REAL getCplus() const {return cplus;}
    REAL getCminus() const {return cminus;}
    REAL getMaxRadius() const;	// get max among aplus, aminus, bplus, bminus, cplus and cminus. August 21, 2013
    REAL getMinRadius() const;	// get min. August 27, 2013
    REAL getYoung() const {return young;}
    REAL getPoisson() const {return poisson;};
    REAL getVolume() const {return volume;}
    REAL getMass() const {return mass;}
    REAL getDensity() const {return density;}
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

    REAL getRadius(vec v, int, int, int) const;	// August 21, 2013
//    REAL getRadius(vec v);
    REAL getTransEnergy() const;
    REAL getRotatEnergy() const;
    REAL getKinetEnergy() const;
    REAL getPotenEnergy(REAL ref) const;	// August 16, 2013
    
    void setID(int n){ID=n;}
    void setType(int n) {type=n;}
    void setAplus(REAL dd){aplus=dd;}	
    void setAminus(REAL dd){aminus=dd;}
    void setBplus(REAL dd){bplus=dd;}	// August 16, 2013
    void setBminus(REAL dd){bminus=dd;}
    void setCplus(REAL dd){cplus=dd;}
    void setCminus(REAL dd){cminus=dd;}
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
    
    void clearForce();	// Add moment by gravity. August 28, 2013
    void addForce(vec vv) {force+=vv;}
    void addMoment(vec vv) {moment+=vv;}
    void update();	// August 16, 2013
    
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
    vec curr_direction_a, curr_direction_b, curr_direction_c;// the direction of the three axles, in radian
    vec prev_direction_a, prev_direction_b, prev_direction_c;
    vec curr_velocity;   // the velocity of the mass center
    vec prev_velocity;
    vec curr_omga;       // angular velocity in global frame!
    vec prev_omga;
    vec curr_acceleration;
    vec prev_acceleration;
    vec pre_force;
    vec force;
    vec pre_moment;
    vec moment;
    vec const_force;
    vec const_moment;
    vec flb_force;
    vec flb_moment;
    vec mres;            // resistence moment provided by normal distribution force
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
    void init();	// August 19, 2013
    
  public:
    int IsFBP;
  };
  
} // namespace dem ends

#endif
