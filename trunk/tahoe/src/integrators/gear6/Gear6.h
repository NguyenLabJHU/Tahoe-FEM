#ifndef _GEAR6_H_
#define _GEAR6_H_

#include "Environment.h"

/* base class */
#include "ControllerT.h"

namespace Tahoe {

/** explicit, central differences time integrator */
class Gear6: virtual public ControllerT
{
public:

	/** constructor */
	Gear6(void) { };

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

#endif /* _GEAR6_H_ */
