/* $Id: TrapezoidIntegrator.h,v 1.2 2002-07-02 19:55:10 cjkimme Exp $ */
/* created: paklein (10/03/1999) */

#ifndef _TRAPEZOID_CONTROLLER_H_
#define _TRAPEZOID_CONTROLLER_H_

/* base classes */
#include "nTrapezoid.h"
#include "eTrapezoid.h"

/* forward declarations */
#include "ios_fwd_decl.h"


namespace Tahoe {

class TrapezoidIntegrator: public nTrapezoid, public eTrapezoid
{
public:

	/* constructor */
	TrapezoidIntegrator(ostream& out);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe 
#endif /* _TRAPEZOID_CONTROLLER_H_ */
