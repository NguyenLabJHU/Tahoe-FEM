/* $Id: D2OrthoMLS2DT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/17/1999)                                          */

#ifndef _D2_ORTHO_MLS_2D_T_H_
#define _D2_ORTHO_MLS_2D_T_H_

/* base class */
#include "D2OrthoMLSSolverT.h"

class D2OrthoMLS2DT: public D2OrthoMLSSolverT
{
public:

	/* constructor */
	D2OrthoMLS2DT(int complete);
	
protected:

	/* return the number of monomial terms for the given completeness */
	virtual int NumberOfMonomials(int completeness) const;

	/* return monomials evaluated at coords */
	virtual void SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp);
	virtual void _SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp,
		dArray2DT& DDp);
};

#endif /* _D2_ORTHO_MLS_2D_T_H_ */
