/* $Id: ParticlePropertyT.h,v 1.1.6.1 2003-04-09 16:01:46 cjkimme Exp $ */
#ifndef _PARTICLE_PROPERTY_T_H_
#define _PARTICLE_PROPERTY_T_H_

#include "ios_fwd_decl.h"

namespace Tahoe {

/** base class for particle properties and interactions */
class ParticlePropertyT
{
public:

	/** enum for particle property types */
	enum TypeT {
        kHarmonicPair = 0, /**< harmonic pair potential */
    kLennardJonesPair = 1, /**< Jennard-Jones 6/12 pair potential */
         kParadynPair = 2, /**< pair potential in Paradyn (EAM) format */
          kParadynEAM = 3  /**< EAM potentials in Paradyn format */
	};
	
	enum ThermostatT {
	    kFreeParticle = 0, /**< you figure it out */
      kDampedParticle = 1, /**< velocity-dependent damping */
    kLangevinParticle = 2  /**< Langevin (stochastic) thermostat */
	};
	
	/** stream extraction operators */
	friend istream& operator>>(istream& in, ParticlePropertyT::TypeT& property);	
	friend istream& operator>>(istream& in, ParticlePropertyT::ThermostatT& property);	

	/** constructor */
	ParticlePropertyT(void);

	/** destructor */
	virtual ~ParticlePropertyT(void) {};
	
	/** interaction distance. Distance used for doing neighbor searches and
	 * determining the depth of interprocessor communication layers. */
	double Range(void) const { return fRange; };

	/** particle mass */
	double Mass(void) const { return fMass; };

	/** write properties to output */
	virtual void Write(ostream& out) const;

protected:

	/** \name methods to set particle properties */
	/*@{*/
	void SetMass(double mass) { fMass = mass; };
	void SetRange(double range) { fRange = range; };
	/*@}*/
	
protected:

	/** \name properties */
	/*@{*/
	double fMass;	
	double fRange;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARTICLE_PROPERTY_T_H_ */
