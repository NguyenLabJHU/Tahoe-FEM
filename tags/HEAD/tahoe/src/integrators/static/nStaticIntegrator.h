/* $Id: nStaticIntegrator.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */

#ifndef _N_STATICCONTROLLER_H_
#define _N_STATICCONTROLLER_H_

/* base classes */
#include "IntegratorT.h"
#include "nIntegratorT.h"

class nStaticIntegrator: public virtual IntegratorT, public nIntegratorT
{
public:

	/* constructor */
	nStaticIntegrator(void);

	/* consistent BC's */
	virtual void ConsistentKBC(const KBC_CardT& KBC);

	/* pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const;
	  	
protected:  	
	
	/* recalculate time stepping constants */
	virtual void nComputeParameters(void);
	
};

#endif /* _N_STATICCONTROLLER_H_ */
