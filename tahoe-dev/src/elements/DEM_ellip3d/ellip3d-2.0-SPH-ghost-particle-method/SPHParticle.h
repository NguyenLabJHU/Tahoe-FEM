#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H


#include <iostream>
#include <vector>
#include <math.h>

#include "vec.h"
#include "realtypes.h"
#include "globfuncs.h"
#include "parameter.h"

namespace dem{
class particle;
}


namespace sph{

class SPHParticle {

public:
    // Default Constructor
    SPHParticle();
    SPHParticle(REAL mass, REAL density, REAL x, REAL y, REAL z, int t);
    SPHParticle(REAL mass, REAL density, REAL x, REAL y, REAL z, dem::vec local, dem::particle* p, int t);

    ~SPHParticle() {};

    void calculateParticleVolume() {volume=mass/density;}
    void calculateParticlePressure() {pressure= dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);}
    void calculateParticleViscosity() {mu=density*dem::nu;}	// dynamic viscosity

    void setDensityDot(REAL a) {densityDot = a;}
    void setVelocityDot(const dem::vec& a) {velocityDot = a;}
    void setDensityDotVelocityDotZero() {densityDot = 0; velocityDot = 0; velocityCorrection = 0;}
    void setCurrPosition(dem::vec a) {curr_x = a;}
    void setCurrPositionX(REAL a) {curr_x.setx(a);}
    void setCurrPositionY(REAL a) {curr_x.sety(a);}
    void setCurrPositionInitial() {curr_x = initial_X;}
    void setCurrVelocity(dem::vec a) {velocity = a;}
    void setType(int a) {type=a;}
    void addVelocityDot(const dem::vec& a) {velocityDot = velocityDot+a;}
    void addVelocityCorrection(const dem::vec& a) {velocityCorrection = velocityCorrection+a;}
    void addDensityDot(REAL a) {densityDot = densityDot+a;}
    void addCurrPositionX(REAL a) {curr_x.setx(curr_x.getx()+a);}

    REAL getParticleMass() const {return mass;}
    REAL getParticleDensity() const {return density;}
    REAL getDensityDot() const {return densityDot;}
    REAL getParticleVolume() const {return volume;}	// everytime when use getParticleVolume(), make sure that volume has been calculated!!!!
    REAL getParticlePressure() const {return pressure;}	// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
    REAL getParticleViscosity() const {return mu;}	// everytime when use getParticleViscosity(), make sure that viscosity has been calculated!!!!
    REAL getKineticEnergy() const {return 0.5*mass*velocity%velocity;}
    dem::vec getInitPosition() const {return initial_X;}
    dem::vec getCurrPosition() const {return curr_x;}
    dem::vec getLocalPosition() const{return local_X;}
    dem::vec getTrialPosition() const {return curr_x+velocity*(dem::TIMESTEP);}	// return the trial position, used for the boundary handling model based on momentum conservation
    dem::vec getDisplacement() const {return curr_x-initial_X;}
    dem::vec getVelocity() const {return velocity;}
    dem::vec getVelocityDot() const {return velocityDot;}
    int getType() const {return type;}
    dem::particle* getDemParticle() {return demParticle;}
	
    void fixYandZ() {velocityDot.sety(0); velocityDot.setz(0); velocityCorrection.sety(0); velocityCorrection.setz(0);}
    void fixZ() {velocityDot.setz(0);velocityCorrection.setz(0);}
    void fixY() {velocityDot.sety(0);velocityCorrection.sety(0);}
    void fixXYZ() {velocityDot=0;velocityCorrection=0;}
    void fixZinApplyBoundary(REAL fix_z) {velocityDot.setz(0); velocityCorrection.setz(0); velocity.setz(0); curr_x.setz(fix_z);}
    void initial();	// initial displacement, velocity and acceleration
    void updateParticle();	// update density, velocity and positions
    void updateParticleDensity();	// only update density

    void updateParticlePositionDensityLeapFrog();
    void updateParticleVelocityLeapFrog();
    void initialParticleVelocityLeapFrog();

	
private:

    REAL mass;		// mass
    REAL density;	// mass density
    REAL volume;	// volume, v=m/density, need to be calculated before being used!!!
    REAL pressure;	// pressure, need to be calculated before being used!!!
    REAL mu;	// dynamic viscosity, mu=density*dem::nu; need to be calculated before being used!!!
    REAL densityDot;	// density dot

    dem::vec curr_x;	// current position
    dem::vec initial_X;	// initial position

    dem::vec velocity;	// velocity of SPH particle
    dem::vec velocityDot;	// velocity dot, acceleration
    dem::vec velocityCorrection;	// velocity correction, the delta_a term, by Monanghan's paper(1994)

    dem::vec local_X;	// for the ghost point only, the local coordinates in the dem particle

    // variable for linked list searching
    int type;	// particle type: 1, free particle; 2, ghost particle; 3, boundary particle
    dem::particle* demParticle;	// if is ghost particle, then pointer to its dem particle
   
}; // end particle


} // end namespace sph

#endif
