/* $Id: D2OrthoMLS1DT.h,v 1.2 2002-07-02 19:57:02 cjkimme Exp $ */
/* created: paklein (10/21/1999)                                          */

#ifndef _D2_ORTHO_MLS_1D_T_H_
#define _D2_ORTHO_MLS_1D_T_H_

/* base class */
#include "D2OrthoMLSSolverT.h"


namespace Tahoe {

class D2OrthoMLS1DT: public D2OrthoMLSSolverT
{
public:

	/* constructor */
	D2OrthoMLS1DT(int complete);
	
protected:

	/* return the number of monomial terms for the given completeness */
	virtual int NumberOfMonomials(int completeness) const;

	/* return monomials evaluated at coords */
	virtual void SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp);
	virtual void _SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp,
		dArray2DT& DDp);
};

} // namespace Tahoe 
#endif /* _D2_ORTHO_MLS_1D_T_H_ */
