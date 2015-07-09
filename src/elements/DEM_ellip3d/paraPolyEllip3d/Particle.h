#ifndef PARTICLE_H
#define PARTICLE_H

#include "Parameter.h"
#include "realtypes.h"
#include "Vec.h"
#include "Gradation.h"
#include "Rectangle.h"
#include "Cylinder.h"
#include "Boundary.h"
#include <cstddef>
#include <map>
#include <vector>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace dem {
  
  class Particle {

  private:
    // types of individual particle:
    //   0 - free particle
    //   1 - fixed particle
    //   2 - special case 2 (pure moment): translate first, then rotate only, MNT_START needs to be defined
    //   3 - special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
    //   4 - special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only
    //   5 - free boundary particle
    //   6 - translate only, no rotation
    //  10 - ghost particle
    std::size_t  id;
    std::size_t  type;            
    REAL aplus, aminus, bplus, bminus, cplus, cminus;   // six parameters, not order by value, corresponding to three principal axels, July 1, 2015.
    REAL young;      // note: a(currDirecA), b(currDirecB), c(currDirecC) corresponds to x, y, z in local frame, respectively
    REAL poisson;
    Vec  currPos;    // particle center
    Vec  prevPos;    // for poly-ellipsoid, the geometry centroid and the mass center is not the same, July 1, 2015
    Vec  localCenterMass;	// local coordinate of mass center with center_geo as origin and princinpal directions as axels, unchanged. July 1, 2015
    Vec  prevCenterMass;
    Vec  currCenterMass;	//global coordinates of previous and current center of mass
    Vec  initCenterMass;	// initial center
    Vec  startCenterMass;	// position at starting step when tessellation is regenerated
    Vec  currDirecA, currDirecB, currDirecC; // direction of the three axles, in radian
    Vec  prevDirecA, prevDirecB, prevDirecC;
    Vec  currVeloc;  // the velocity of the mass center
    Vec  prevVeloc;
    Vec  currOmga;   // angular velocity in global frame!
    Vec  prevOmga;
    Vec  force;
    Vec  prevForce;
    Vec  moment;
    Vec  prevMoment;
    Vec  constForce;
    Vec  constMoment;
    REAL density;    // specific gravity
    REAL mass;
    REAL volume;
    Vec  momentJ;    // moment of inertia in local body-fixed frame
    REAL coef1[10];// record particle's coefficients in global coordinates
    REAL coef2[10];
    REAL coef3[10];
    REAL coef4[10];
    REAL coef5[10];	// July 1, 2015
    REAL coef6[10];
    REAL coef7[10];
    REAL coef8[10];
    REAL kinetEnergy;// kinetic energy
    std::size_t  contactNum;
    bool inContact;  // in contact with other particle or boundary
    std::vector< std::vector<REAL> > fluidGrid;

  public:
    Particle();
    Particle(std::size_t n, std::size_t type, Vec center, REAL r, REAL young, REAL poisson);	// July 3, 2015
    Particle(std::size_t n, std::size_t type, Vec center, REAL aplus, REAL aminus, REAL bplus, REAL bminus, REAL cplus, REAL cminus, REAL young, REAL poisson);	// July 3, 2015
    Particle(std::size_t n, std::size_t type, Vec center, Gradation& grad, REAL young, REAL poisson);	// July 3, 2015
    Particle(std::size_t n, std::size_t type, REAL aplus, REAL aminus, REAL bplus, REAL bminus, REAL cplus, REAL cminus, Vec position, Vec dirca, Vec dircb, Vec dircc, REAL young, REAL poisson);	// July 3, 2015
    
    std::size_t  getId() const {return id;}
    std::size_t  getType() const {return type;}
    REAL getAplus() const {return aplus;}
    REAL getAminus() const {return aminus;}
    REAL getBplus() const {return bplus;}
    REAL getBminus() const {return bminus;}
    REAL getCplus() const {return cplus;}
    REAL getCminus() const {return cminus;}
    REAL getA() const {return 0.5*(aplus+aminus);}	// July 6, 2015. For fluid.cpp
    REAL getB() const {return 0.5*(bplus+bminus);}
    REAL getC() const {return 0.5*(cplus+cminus);}
    REAL getMaxRadius() const;	// get max among aplus, aminus, bplus, bminus, cplus and cminus. July 1, 2015
    REAL getMinRadius() const;	// get min. July 1, 2015
    REAL getYoung() const {return young;}
    REAL getPoisson() const {return poisson;};
    REAL getVolume() const {return volume;}
    REAL getMass() const {return mass;}
    REAL getDensity() const {return density;}
    Vec  getCurrPos() const {return currPos;}
    Vec  getPrevPos() const {return prevPos;}
    Vec  getInitCenterMass() const {return initCenterMass;}	// initial position for granular strain
    Vec  getStartCenterMass() const {return startCenterMass;}
    Vec  getCurrCenterMass() const {return currCenterMass;}	// July 1, 2015
    Vec  getPrevCenterMass() const {return prevCenterMass;}
    Vec  getCurrDirecA() const {return currDirecA;}
    Vec  getCurrDirecB() const {return currDirecB;}
    Vec  getCurrDirecC() const {return currDirecC;}
    Vec  getPrevDirecA() const {return prevDirecA;}
    Vec  getPrevDirecB() const {return prevDirecB;}
    Vec  getPrevDirecC() const {return prevDirecC;}
    Vec  getCurrVeloc() const {return currVeloc;}
    Vec  getPrevVeloc() const {return prevVeloc;}
    Vec  getCurrOmga() const {return currOmga;}
    Vec  getPrevOmga() const {return prevOmga;}
    Vec  getForce() const {return force;}
    Vec  getMoment() const {return moment;}
    Vec  getAccel() const {return force/mass;}
    Vec  getConstForce() const {return constForce;}
    Vec  getConstMoment() const {return constMoment;}
    Vec  getmomentJ() const {return momentJ;}
    bool isInContact() const {return inContact;}
    std::size_t  getContactNum() const {return contactNum;}

    REAL getRadius(Vec v, int, int, int) const;	// July 1, 2015
    REAL getTransEnergy() const;
    REAL getRotatEnergy() const;
    REAL getKinetEnergy() const;
    REAL getPotenEnergy(REAL ref) const;	// July 1, 2015
    
    void setId(std::size_t n) {id = n;}
    void setType(std::size_t n) {type = n;}
    void setAplus(REAL dd){aplus=dd;}	
    void setAminus(REAL dd){aminus=dd;}
    void setBplus(REAL dd){bplus=dd;}	// July 1, 2015
    void setBminus(REAL dd){bminus=dd;}
    void setCplus(REAL dd){cplus=dd;}
    void setCminus(REAL dd){cminus=dd;}
    void expand(REAL percent) {aplus *= (1+percent); aminus *= (1+percent); bplus *= (1+percent); bminus *= (1+percent); cplus *= (1+percent); cminus *= (1+percent);}	// July 1, 2015, actually, the geometry information, such as mass, volume moment inertia should be recalculated after expansion
    void setCurrPos(Vec vv) {currPos = vv;}
    void setPrevPos(Vec vv) {prevPos = vv;}
    void setCurrCenterMass(Vec vv){currCenterMass=vv;}	// July 1, 2015
    void setPrevCenterMass(Vec vv){prevCenterMass=vv;}
    void setInitCenterMass(Vec vv){initCenterMass=vv;}	// initial center for granular strain
    void setStartCenterMass(Vec vv){startCenterMass=vv;}	// position at starting step when tessellation is regenerated
    void setCurrDirecA(Vec vv) {currDirecA = vv;}
    void setCurrDirecB(Vec vv) {currDirecB = vv;}
    void setCurrDirecC(Vec vv) {currDirecC = vv;}
    void setPrevDirecA(Vec vv) {prevDirecA = vv;}
    void setPrevDirecB(Vec vv) {prevDirecB = vv;}
    void setPrevDirecC(Vec vv) {prevDirecC = vv;}
    void setCurrVeloc(Vec vv) {currVeloc = vv;}
    void setPrevVeloc(Vec vv) {prevVeloc = vv;}
    void setCurrOmga(Vec vv) {currOmga = vv;}
    void setPrevOmga(Vec vv) {prevOmga = vv;}
    void setForce(Vec vv) {force = vv;}
    void setMoment(Vec vv) {moment = vv;}
    void setConstForce(Vec vv) {constForce = vv;}
    void setConstMoment(Vec vv) {constMoment = vv;}
    void setmomentJ(Vec v) {momentJ = v;}
    void setMass(REAL d) {mass = d;}
    void setDensity(REAL dn) {density = dn;}
    void setInContact(bool value) {inContact = value;}
    void setContactNum(std::size_t num) {contactNum = num;}

    void clearContactForce();
    void addForce(Vec vv) {force += vv;}
    void addMoment(Vec vv) {moment += vv;}
    void update();	// July 1, 2015
    void dragForce();

    Vec globalToLocal(Vec input) const;
    Vec localToGlobal(Vec input) const;
    
    // update global coefficients in the following form based on position/dimensions/orientations
    // a0 x^2 + a1 y^2 + a2 z^2 + a3 xy + a4 yz + a5 zx + a6 x + a7 y + a8 z + a9 = 0.. July 1, 2015
    void globalCoef();  
    void getGlobalCoef(REAL coef[], int num_oct) const; // retrieve global coeffs into coef[]
    REAL surfaceError(Vec pt) const;	// related to coef, but only used in Fluid.cpp. need to modify it later for poly-ellipsoids. July 1, 2015
    
    // v is the point the line passing through, dirc is the unit vector parallel to the line
    bool intersectWithLine(Vec v, Vec dirc, Vec rt[], int num_oct) const;	// July 1, 2015
    
    // find the point on plane which is deepest into a particles, px + qy + rz + s = 0 is the equation 
    // of the plane, true means intersection; false means no intersection.
    bool nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, Vec &ptnp, int num_oct) const;	// July 1, 2015
    
    // calculate the normal force between particle and a plane rigid boundary
    void planeRBForce(planeBoundary *plane,	// apply moment on the mass center. July 1, 2015
		      std::map<std::size_t,std::vector<BoundaryTgt> > &BoundarytgtMap,
		      std::vector<BoundaryTgt> &vtmp);
    
    // calculate the normal force between particle and a cylinder wall
    Vec cylinderRBForce(std::size_t boundaryId, const Cylinder &S, int side);	// July 3, 2015
    void clearFluidGrid();
    void recordFluidGrid(std::size_t i, std::size_t j, std::size_t k);
    std::vector< std::vector<REAL> > & getFluidGrid() { return fluidGrid; }
    
  private:
    void init();    

  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & id;
      ar & type;
      ar & aplus & aminus & bplus & bminus & cplus & cminus;
      ar & young;
      ar & poisson;
      ar & currPos;
      ar & prevPos;
      ar & localCenterMass;	
      ar & prevCenterMass;
      ar & currCenterMass;	
      ar & initCenterMass;	
      ar & startCenterMass;
      ar & currDirecA & currDirecB & currDirecC;
      ar & prevDirecA & prevDirecB & prevDirecC;
      ar & currVeloc;
      ar & prevVeloc;
      ar & currOmga;
      ar & prevOmga;
      ar & force;
      ar & prevForce;
      ar & moment;
      ar & prevMoment;
      ar & constForce;
      ar & constMoment;
      ar & density;
      ar & mass;
      ar & volume;
      ar & momentJ;
      ar & coef1;
      ar & coef2;
      ar & coef3;
      ar & coef4;
      ar & coef5;
      ar & coef6;
      ar & coef7;
      ar & coef8;
      ar & kinetEnergy;
      ar & contactNum;
      ar & inContact;
      ar & fluidGrid;
    }  
    
  };
  
} // namespace dem ends

#endif
