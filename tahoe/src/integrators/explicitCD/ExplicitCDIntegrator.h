/* $Id: ExplicitCDIntegrator.h,v 1.4.56.1 2004-04-08 07:33:40 paklein Exp $ */
/* created: paklein (03/23/1997) */
#ifndef _EXP_CD_CONTROLLER_H_
#define _EXP_CD_CONTROLLER_H_

#include "Environment.h"

/* base classes */
#include "nExplicitCD.h"
#include "eExplicitCD.h"

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/** controller for an explicit 2nd order accurate, central difference
 * time-stepping algorithm */
class ExplicitCDIntegrator: public nExplicitCD, public eExplicitCD
{
public:

	/** constructor */
	ExplicitCDIntegrator(void);
	  	
protected:  	
	
	/** recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe 
#endif /* _EXP_CD_CONTROLLER_H_ */
