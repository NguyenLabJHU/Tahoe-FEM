/* $Id: ExplicitCDIntegrator.h,v 1.2 2002-04-02 23:19:23 paklein Exp $ */
/* created: paklein (03/23/1997) */

#ifndef _EXP_CD_CONTROLLER_H_
#define _EXP_CD_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nExplicitCD.h"
#include "eExplicitCD.h"

/* forward declarations */
#include "ios_fwd_decl.h"

/** controller for an explicit 2nd order accurate, central difference
 * time-stepping algorithm */
class ExplicitCDIntegrator: public nExplicitCD, public eExplicitCD
{
public:

	/** constructor */
	ExplicitCDIntegrator(ostream& out);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

#endif /* _EXP_CD_CONTROLLER_H_ */
