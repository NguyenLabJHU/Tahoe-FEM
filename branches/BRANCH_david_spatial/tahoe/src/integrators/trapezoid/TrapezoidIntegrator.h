/* $Id: TrapezoidIntegrator.h,v 1.4 2004-07-15 08:30:53 paklein Exp $ */
/* created: paklein (10/03/1999) */
#ifndef _TRAPEZOID_CONTROLLER_H_
#define _TRAPEZOID_CONTROLLER_H_

/* base classes */
#include "nTrapezoid.h"
#include "eTrapezoid.h"

namespace Tahoe {

class TrapezoidIntegrator: public nTrapezoid, public eTrapezoid
{
public:

	/* constructor */
	TrapezoidIntegrator(void);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe 
#endif /* _TRAPEZOID_CONTROLLER_H_ */
