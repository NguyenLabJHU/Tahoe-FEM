/* $Id: FDCubicT.h,v 1.2 2001-07-03 01:35:06 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _FD_CUBIC_T_H_
#define _FD_CUBIC_T_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "CubicT.h"

class FDCubicT: public FDHookeanMatT, public CubicT
{
public:

	/* constructor */
	FDCubicT(ifstreamT& in, const FiniteStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

#endif /* _FD_CUBIC_T_H_ */
