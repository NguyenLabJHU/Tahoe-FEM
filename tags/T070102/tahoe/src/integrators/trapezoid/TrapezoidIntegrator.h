/* $Id: TrapezoidIntegrator.h,v 1.1 2002-04-02 23:20:42 paklein Exp $ */
/* created: paklein (10/03/1999) */

#ifndef _TRAPEZOID_CONTROLLER_H_
#define _TRAPEZOID_CONTROLLER_H_

/* base classes */
#include "nTrapezoid.h"
#include "eTrapezoid.h"

/* forward declarations */
#include "ios_fwd_decl.h"

class TrapezoidIntegrator: public nTrapezoid, public eTrapezoid
{
public:

	/* constructor */
	TrapezoidIntegrator(ostream& out);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

#endif /* _TRAPEZOID_CONTROLLER_H_ */
