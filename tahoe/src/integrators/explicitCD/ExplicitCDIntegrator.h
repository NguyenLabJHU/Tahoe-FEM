/* $Id: ExplicitCDIntegrator.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Integrator for an explicit 2nd order accurate, central difference      */
/* time-stepping algorithm.                                               */

#ifndef _EXP_CD_CONTROLLER_H_
#define _EXP_CD_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nExplicitCD.h"
#include "eExplicitCD.h"

/* forward declarations */
#include "ios_fwd_decl.h"

class ExplicitCDIntegrator: public nExplicitCD, public eExplicitCD
{
public:

	/* constructor */
	ExplicitCDIntegrator(ostream& out);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

#endif /* _EXP_CD_CONTROLLER_H_ */
