/* $Id: ParticlePropertyT.h,v 1.9 2004-07-15 08:29:49 paklein Exp $ */
#ifndef _PARTICLE_PROPERTY_T_H_
#define _PARTICLE_PROPERTY_T_H_

/* base class */
#include "ParameterInterfaceT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/** base class for particle properties and interactions */
class ParticlePropertyT: public ParameterInterfaceT
{
public:

	/** enum for particle property types */
	enum TypeT {
        kHarmonicPair = 0, /**< harmonic pair potential */
    kLennardJonesPair = 1, /**< Jennard-Jones 6/12 pair potential */
         kParadynPair = 2, /**< pair potential in Paradyn (EAM) format */
          kParadynEAM = 3, /**< EAM potentials in Paradyn format */
	  kMatsuiPair = 4  /**< Matsui pair potential */
	};

	/** stream extraction operators */
//	friend istream& operator>>(istream& in, ParticlePropertyT::TypeT& property);

	/** constructor */
	ParticlePropertyT(void);

	/** destructor */
	virtual ~ParticlePropertyT(void) {};
	
	/** interaction distance. Distance used for doing neighbor searches and
	 * determining the depth of interprocessor communication layers. */
	double Range(void) const { return fRange; };

	/** particle mass */
	double Mass(void) const { return fMass; };

	/** nominal nearest neighbor distance */
	double NearestNeighbor(void) const { return fNearestNeighbor; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** \name methods to set particle properties */
	/*@{*/
	void SetMass(double mass) { fMass = mass; };
	void SetRange(double range) { fRange = range; };
	void SetNearestNeighbor(double nearest) { fNearestNeighbor = nearest; };
	/*@}*/
	
protected:

	/** \name properties */
	/*@{*/
	double fMass;	
	double fRange;
	double fNearestNeighbor;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PARTICLE_PROPERTY_T_H_ */
