/* $Id: StaticIntegrator.h,v 1.5 2002-07-05 22:27:55 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _STATIC_CONTROLLER_H_
#define _STATIC_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nStaticIntegrator.h"
#include "eStaticIntegrator.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

class StaticIntegrator: public nStaticIntegrator, public eStaticIntegrator
{
public:

	/** constructor */
	StaticIntegrator(ostream& out);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
	
};

} // namespace Tahoe 
#endif /* _STATIC_CONTROLLER_H_ */
