#ifndef _VERLET_CONTROLLER_H_
#define _VERLET_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nVerlet.h"
#include "eVerlet.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/** controller for an explicit 4th order accurate velocity verlet
 * time integration algorithm */
class VerletIntegrator: public nVerlet, public eVerlet
{
public:

	/** constructor */
	VerletIntegrator(ostream& out);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe

#endif /* _VERLET_CONTROLLER_H_ */
