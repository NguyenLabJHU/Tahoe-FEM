/* $Id: ExplicitCD.h,v 1.2 2002-07-02 19:55:08 cjkimme Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _EXPLICIT_CD_H_
#define _EXPLICIT_CD_H_

#include "Environment.h"

/* base class */
#include "ControllerT.h"

/** explicit, central differences time integrator */

namespace Tahoe {

class ExplicitCD: virtual public ControllerT
{
public:

	/** constructor */
	ExplicitCD(void) { };

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
#endif /* _EXPLICIT_CD_H_ */
