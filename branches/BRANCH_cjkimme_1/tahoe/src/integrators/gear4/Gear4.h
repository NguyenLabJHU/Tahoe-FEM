#ifndef _GEAR_04_H_
#define _GEAR_04_H_

#include "Environment.h"

/* base class */
#include "IntegratorT.h"

namespace Tahoe {

/** Gear time integrator. Base class for the Gear predictor-corrector 
 * time integrator with a local trunction error of 
 * \f$ \Delta t^4 \f$. This is implemented to work with barostats
 * for ParticleT elements. */
class Gear4: virtual public IntegratorT
{
public:

	/** constructor */
	Gear4(void) { };

	/** \name integrator parameters */
	/*@{*/
	/** return flag indicating whether integrator is implicit or explicit */
	virtual ImpExpFlagT ImplicitExplicit(void) const { return kExplicit; };

	/** return order time discretization */
	virtual int Order(void) const { return 3; };

	/** return order field derivative which is treated as the primary 
	 * unknown value */
	virtual int OrderOfUnknown(void) const { return 2; };
	/*@}*/
};

} // namespace Tahoe

#endif /* _GEAR_04_H_ */
