#ifndef _GEAR_04_INTEGRATOR_H_
#define _GEAR_06_INTEGRATOR_H_

#include "Environment.h"

/* base classes */
#include "nGear4.h"
#include "eGear4.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/** controller for an explicit 6th order Gear predictor corrector
 * time integration algorithm */
class Gear4Integrator: public nGear4, public eGear4
{
public:

	/** constructor */
	Gear4Integrator(ostream& out);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe

#endif /* _GEAR_04_CONTROLLER_H_ */
