/* $Id: StaticIntegrator.h,v 1.6 2004-07-15 08:30:48 paklein Exp $ */
/* created: paklein (10/14/1996) */
#ifndef _STATIC_CONTROLLER_H_
#define _STATIC_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nStaticIntegrator.h"
#include "eStaticIntegrator.h"

namespace Tahoe {

class StaticIntegrator: public nStaticIntegrator, public eStaticIntegrator
{
public:

	/** constructor */
	StaticIntegrator(void);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
	
};

} // namespace Tahoe 
#endif /* _STATIC_CONTROLLER_H_ */
