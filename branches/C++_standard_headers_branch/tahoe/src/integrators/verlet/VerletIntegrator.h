/* $Id: VerletIntegrator.h,v 1.3 2004-07-15 08:30:57 paklein Exp $ */
#ifndef _VERLET_CONTROLLER_H_
#define _VERLET_CONTROLLER_H_

/* base classes */
#include "nVerlet.h"
#include "eVerlet.h"

namespace Tahoe {

/** controller for an explicit 4th order accurate velocity verlet
 * time integration algorithm */
class VerletIntegrator: public nVerlet, public eVerlet
{
public:

	/** constructor */
	VerletIntegrator(void);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe

#endif /* _VERLET_CONTROLLER_H_ */
