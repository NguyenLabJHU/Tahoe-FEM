/* $Id: StaticT.h,v 1.3 2002-07-05 22:27:55 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _STATIC_T_H_
#define _STATIC_T_H_

/* base class */
#include "ControllerT.h"

namespace Tahoe {

/** explicit, central differences time integrator */
class StaticT: public virtual ControllerT
{
public:

	/** constructor */
	StaticT(void) { };

	/** \name integrator parameters */
	/*@{*/
	/** return flag indicating whether integrator is implicit or explicit */
	virtual ImpExpFlagT ImplicitExplicit(void) const { return kImplicit; };

	/** return order time discretization */
	virtual int Order(void) const { return 0; };

	/** return order field derivative which is treated as the primary 
	 * unknown value */
	virtual int OrderOfUnknown(void) const { return 0; };
	/*@}*/
};

} // namespace Tahoe 
#endif /* _STATIC_T_H_ */
