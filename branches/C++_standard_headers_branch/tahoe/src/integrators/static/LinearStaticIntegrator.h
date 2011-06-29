/* $Id: LinearStaticIntegrator.h,v 1.3 2004-07-15 08:30:48 paklein Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _LINEAR_STATIC_CONTROLLER_H_
#define _LINEAR_STATIC_CONTROLLER_H_

/* base classes */
#include "nLinearStaticIntegrator.h"
#include "eStaticIntegrator.h"

namespace Tahoe {

class LinearStaticIntegrator: public nLinearStaticIntegrator, public eStaticIntegrator
{
public:

	/** constructor */
	LinearStaticIntegrator(void);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
	
};

} // namespace Tahoe 
#endif /* _LINEAR_STATIC_CONTROLLER_H_ */
