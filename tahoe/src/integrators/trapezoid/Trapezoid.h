/* $Id: Trapezoid.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/03/1999)                                          */

#ifndef _TRAPEZOID_H_
#define _TRAPEZOID_H_

/* base classes */
#include "nTrapezoid.h"
#include "eTrapezoid.h"

/* forward declarations */
#include "ios_fwd_decl.h"

class Trapezoid: public nTrapezoid, public eTrapezoid
{
public:

	/* constructor */
	Trapezoid(ostream& out);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

#endif /* _TRAPEZOID_H_ */
