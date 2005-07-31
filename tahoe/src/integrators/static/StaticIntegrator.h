/* $Id: StaticIntegrator.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */
/* This controller does not have a node controller branch                 */

#ifndef _STATICCONTROLLER_H_
#define _STATICCONTROLLER_H_

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

#endif /* _STATICCONTROLLER_H_ */
