/* $Id: Verlet.h,v 1.3.46.1 2004-04-08 07:33:46 paklein Exp $ */
#ifndef _VERLET_H_
#define _VERLET_H_

/* base class */
#include "IntegratorT.h"

namespace Tahoe {

/** explicit, central differences time integrator */
class Verlet: virtual public IntegratorT
{
public:

	/** constructor */
	Verlet(void) { };

	/** \name integrator parameters */
	/*@{*/
	/** return flag indicating whether integrator is implicit or explicit */
	virtual ImpExpFlagT ImplicitExplicit(void) const { return kExplicit; };

	/** return order time discretization */
	virtual int Order(void) const { return 2; };

	/** return order field derivative which is treated as the primary 
	 * unknown value */
	virtual int OrderOfUnknown(void) const { return 2; };
	/*@}*/
};

} // namespace Tahoe

#endif /* _VERLET_H_ */
