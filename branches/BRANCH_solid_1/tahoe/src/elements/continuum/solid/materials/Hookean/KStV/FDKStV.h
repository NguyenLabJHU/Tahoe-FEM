/* $Id: FDKStV.h,v 1.1.1.1.2.2 2001-06-22 14:18:04 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _FD_KSTV_H_
#define _FD_KSTV_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "IsotropicT.h"

class FDKStV: public FDHookeanMatT, public IsotropicT
{
public:

	/* constructor */
	FDKStV(ifstreamT& in, const FiniteStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

#endif /* _FD_KSTV_H_ */
