/* $Id: LinearStaticIntegrator.h,v 1.1 2001-08-27 17:12:14 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _LINEAR_STATIC_CONTROLLER_H_
#define _LINEAR_STATIC_CONTROLLER_H_

/* base classes */
#include "nLinearStaticIntegrator.h"
#include "eStaticIntegrator.h"

class LinearStaticIntegrator: public nLinearStaticIntegrator, public eStaticIntegrator
{
public:

	/** constructor */
	LinearStaticIntegrator(ostream& out);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
	
};

#endif /* _LINEAR_STATIC_CONTROLLER_H_ */
