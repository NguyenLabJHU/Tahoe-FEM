/* $Id: Trapezoid.h,v 1.2 2002-04-02 23:19:25 paklein Exp $ */

#ifndef _TRAPEZOID_H_
#define _TRAPEZOID_H_

/* base class */
#include "ControllerT.h"

/** implicit, first-order time integrator */
class Trapezoid: virtual public ControllerT
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

#endif /* _TRAPEZOID_H_ */
