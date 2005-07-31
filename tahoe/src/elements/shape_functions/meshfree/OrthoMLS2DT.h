/* $Id: OrthoMLS2DT.h,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (07/03/1998)                                          */

#ifndef _ORTHO_MLS_2D_T_H_
#define _ORTHO_MLS_2D_T_H_

/* base class */
#include "OrthoMLSSolverT.h"

class OrthoMLS2DT: public OrthoMLSSolverT
{
public:

	/* constructor */
	OrthoMLS2DT(int complete);
	
protected:

	/* return the number of monomial terms for the given completeness */
	virtual int NumberOfMonomials(int completeness) const;

	/* evaluate monomials and derivatives at coords */
	virtual void SetMonomials(const dArrayT& coords, dArrayT& p,  dArray2DT& Dp);

};

#endif /* _ORTHO_MLS_2D_T_H_ */
