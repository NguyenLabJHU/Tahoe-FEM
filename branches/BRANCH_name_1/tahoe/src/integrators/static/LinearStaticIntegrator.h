/* $Id: LinearStaticIntegrator.h,v 1.1.6.1 2002-06-27 18:02:31 cjkimme Exp $ */
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
	LinearStaticIntegrator(ostream& out);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
	
};

} // namespace Tahoe 
#endif /* _LINEAR_STATIC_CONTROLLER_H_ */
