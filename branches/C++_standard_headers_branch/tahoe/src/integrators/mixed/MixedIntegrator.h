/* $Header: /home/cvs/t/ta/tahoe/tahoe/src/integrators/mixed/MixedIntegrator.h,v 1.2 2006-08-18 01:15:50 a-kopacz Exp $ */
/* created: a-kopacz (08/08/2006) */

#ifndef _MIXED_CONTROLLER_H_
#define _MIXED_CONTROLLER_H_

/* base classes */
#include "nMixed.h"
#include "eMixed.h"

namespace Tahoe {

class MixedIntegrator: public nMixed, public eMixed
{
public:

	/* constructor */
	MixedIntegrator(void);
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void ComputeParameters(void);
};

} // namespace Tahoe 
#endif /* _MIXED_CONTROLLER_H_ */
