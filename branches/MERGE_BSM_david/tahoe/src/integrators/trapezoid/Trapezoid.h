/* $Id: Trapezoid.h,v 1.6 2004-07-15 08:30:53 paklein Exp $ */
#ifndef _TRAPEZOID_H_
#define _TRAPEZOID_H_

/* base class */
#include "IntegratorT.h"

namespace Tahoe {

/** implicit, first-order time integrator */
class Trapezoid: virtual public IntegratorT
{
public:

	/** constructor */
	Trapezoid(void) { };

	/** \name integrator parameters */
	/*@{*/
	/** return flag indicating whether integrator is implicit or explicit */
	virtual ImpExpFlagT ImplicitExplicit(void) const { return kImplicit; };

	/** return order time discretization */
	virtual int Order(void) const { return 1; };

	/** return order field derivative which is treated as the primary 
	 * unknown value */
	virtual int OrderOfUnknown(void) const { return 1; };
	/*@}*/
};

} // namespace Tahoe 
#endif /* _TRAPEZOID_H_ */
