/* $Id: nIntegratorT.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (10/14/1996)                                          */
/* General interface for a controller.  This is the interface between the NodeManagerT */
/* and the Controller class heirarchy.                                    */

#ifndef _N_CONTROLLERT_H_
#define _N_CONTROLLERT_H_

/* needed for CodeT enum */
#include "KBC_CardT.h"

/* forward declarations */
class dArray2DT;

class nIntegratorT
{
public:

	/* constructor */
	nIntegratorT(void);

	/* destructor */
	virtual ~nIntegratorT(void);

	/* register field arrays */
	virtual void SetField(dArray2DT& field, int order);

	/* consistent BC's */
	virtual void ConsistentKBC(const KBC_CardT& KBC) = 0;

	/* pseudo-boundary conditions for external nodes */
	virtual KBC_CardT::CodeT ExternalNodeCondition(void) const = 0;

protected:  	
	
	/* recalculate time stepping constants */
	virtual void nComputeParameters(void) = 0;
	
protected:

	/* 0th order field data */
	dArray2DT* fU;
};

#endif /* _N_CONTROLLERT_H_ */
