/* $Id: CSEIsoT.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (11/19/1997)                                          */
/* Cohesive surface elements with scalar traction potentials,             */
/* i.e., the traction potential is a function of the gap magnitude,       */
/* or effective gap magnitude only.                                       */

#ifndef _CSE_ISO_T_H_
#define _CSE_ISO_T_H_

/* base class */
#include "CSEBaseT.h"

/* direct members */
#include "pArrayT.h"

/* forward declarations */
class C1FunctionT;

class CSEIsoT: public CSEBaseT
{
public:

	/* constructor */
	CSEIsoT(FEManagerT& fe_manager);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* initialize class data */
	virtual void Initialize(void);

protected:

	/* tangent matrix */
	virtual void LHSDriver(void);

	/* force vector */
	virtual void RHSDriver(void);

	/* driver for nodal value calculations */
	virtual void ComputeNodalValues(const iArrayT& codes);
	
protected:

	/* cohesive surface potentials */
	pArrayT<C1FunctionT*> fSurfPots;
};

#endif /* _CSE_ISO_T_H_ */
