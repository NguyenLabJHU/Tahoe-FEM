/* $Id: ParticlePropertyT.h,v 1.1 2002-11-25 07:19:46 paklein Exp $ */
#ifndef _PARTICLE_PROPERTY_T_H_
#define _PARTICLE_PROPERTY_T_H_

#include "ios_fwd_decl.h"

namespace Tahoe {

/** base class for particle properties and interactions */
class ParticlePropertyT
{
public:

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
