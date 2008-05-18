#ifndef PARTICLE_H
#define PARTICLE_H

#include <map>
#include "vec.h"
#include "gradation.h"
#include "contact.h"
#include "boundary.h"
#include "rectangle.h"
#include "cylinder.h"
#include "boundarytgt.h"

namespace dem {

class particle{
	friend class assembly;
	friend class contact<particle>;
	friend class flb_bdry<particle>;
public:
	particle(int n, int type, vec center, long double r);
	particle(int n, int type, vec center, long double a, long double b, long double c);
	particle(int id,int type, vec dim, vec position, vec dirca, vec dircb, vec dircc);
	particle(int n, int type, vec center, gradation& grad);

	int    getID() const;
	int    getType() const;
	long double getA() const;
	long double getB() const;
	long double getC() const;
	long double getRadius(vec v) const;
	long double getVolume() const;
	long double getMass() const;
	long double getDensity() const;
	vec    getJ() const;
	vec    getCurrPosition() const;
	vec    getPrevPosition() const;
	vec    getCurrDirecA() const;
	vec    getCurrDirecB() const;
	vec    getCurrDirecC() const;
	vec    getPrevDirecA() const;
	vec    getPrevDirecB() const;
	vec    getPrevDirecC() const;
	vec    getCurrVelocity() const;
	vec    getPrevVelocity() const;
	vec    getCurrOmga() const;
	vec    getPrevOmga() const;
	vec    getCurrAcceleration() const;
	vec    getPrevAcceleration() const;
	vec    getCurrAlf() const;
	vec    getPrevAlf() const;
	vec    getForce() const;
	vec    getMoment() const;
	vec    getConstForce() const;
	vec    getConstMoment() const;
	long double getTransEnergy() const;
	long double getRotaEnergy() const;
	long double getKinEnergy() const;
	long double getPotEnergy(long double ref) const;

	void   setID(int n);
	void   setType(int n);
	void   setA(long double dd);
	void   setB(long double dd);
	void   setC(long double dd);
	void   setMass(long double d);
	void   setJ(vec vv);
	void   setCurrPosition(vec vv);
	void   setPrevPosition(vec vv);
	void   setCurrDirecA(vec vv);
	void   setCurrDirecB(vec vv);
	void   setCurrDirecC(vec vv);
	void   setPrevDirecA(vec vv);
	void   setPrevDirecB(vec vv);
	void   setPrevDirecC(vec vv);
	void   setCurrVelocity(vec vv);
	void   setPrevVelocity(vec vv);
	void   setCurrOmga(vec vv);
	void   setPrevOmga(vec vv);
	void   setCurrAcceleration(vec vv);
	void   setPrevAcceleration(vec vv);
	void   setCurrAlf(vec vv);
	void   setPrevAlf(vec vv);
	void   setDensity(long double dn) {density=dn;}
	void   setExternForce(vec fc) {force=fc;}
	void   setExternMoment(vec mm) {moment=mm;}
	void   setForce(vec vv);
	void   setMoment(vec vv);
	void   setConstForce(vec vv);
	void   setConstMoment(vec vv);

	void   setZero();
	void   addForce(vec vv) {force+=vv;}
	void   addMoment(vec vv) {moment+=vv;}
	void   update();

	// update global coefficients in the following form based on position/dimensions/orientations
	// a0 x^2 + a1 y^2 + a2 z^2 + a3 xy + a4 yz + a5 zx + a6 x + a7 y + a8 z + a9 = 0
	void   GlobCoef();  
	void   getGlobCoef(long double coef[]) const; // fetch global coeffs into coef[]
	long double surfaceError(vec pt) const;
	vec    localVec(vec) const;   // transform a vector in global coordinates into local coordinates
	vec    globalVec(vec) const;  // transform a vector in local coordinates into global coordinates
	void   print() const;

	// Assumption: a particle only intersects a plane a little and it cannot pass through the plane
        //             with most of its body, this is guaranteed by contacting forces.

	//v is the point the line passing through
	//d is the unit vector parallel to the line
	bool intersectWithLine(vec v, vec dirc, vec rt[]) const;

	//find the point on plane which is deepest into a particles, px+qy+rz+s=0 is the equation of the plane
	//true means intersection; false means no intersection.
	bool nearestPTOnPlane(long double p, long double q, long double r, long double s, vec& ptnp) const;

	//side indicates which side the particles are in about the plane
	//calculate the normal force between particle and a plane rigid boundary
	void planeRBForce(plnrgd_bdry<particle>* plb,
			  std::map<int,std::vector<boundarytgt> >& BoundarytgtMap,
			  std::vector<boundarytgt>& vtmp,
			  long double &penetr);

	//if side<0, particles are inside the cylinder; else side is outside the cylinder
	//calculate the normal force between particle and a cylinder wall
	vec cylinderRBForce(int bdry_id, const cylinder& S, int side);

private:
        // particle type:
	// 0-free
        // 1-fixed
        // 2-special case 1 (pure moment): translate first, then rotate only, MNT_START needs to be defined
        // 3-special case 2 (ellipsoidal pile): translation in vertical direction only
	int ID;
	int type;            
	long double a,b,c;   // three semi-axles, a>=b>=c
	vec curr_position;   // particle center
	vec prev_position;
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
	long double density; // specific gravity
	long double mass;
	long double volume;
	vec  J;              // moment of inertia in local body-fixed frame
	long double coef[10];// records particle's coefficients in global coordinates
	long double kinEnergy; //kinetic energy
	int  cntnum;

public:
	int IsFBP;
};

} // namespace dem ends

#endif
