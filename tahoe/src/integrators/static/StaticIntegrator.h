/* $Id: StaticIntegrator.h,v 1.2 2001-08-27 17:12:15 paklein Exp $ */
/* created: paklein (10/14/1996) */

#ifndef _STATIC_CONTROLLER_H_
#define _STATIC_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nStaticIntegrator.h"
#include "eStaticIntegrator.h"

/* forward declarations */
#include "ios_fwd_decl.h"

class StaticIntegrator: public nStaticIntegrator, public eStaticIntegrator
{
public:

	/* constructor */
	StaticIntegrator(ostream& out);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
	
};

#endif /* _STATIC_CONTROLLER_H_ */
