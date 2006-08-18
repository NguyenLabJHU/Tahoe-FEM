/* $Header: */
/* created: a-kopacz (08/08/2006) */

#ifndef _MIXED_H_
#define _MIXED_H_

/* base class */
#include "IntegratorT.h"

namespace Tahoe {

/** implicit, first-order time integrator */
class Mixed: virtual public IntegratorT
{
public:

	/** constructor */
	Mixed(void) { };

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
#endif /* _MIXED_H_ */
