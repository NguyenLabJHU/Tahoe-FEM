/* $Id: FDKStV.h,v 1.3.8.1 2002-10-28 06:48:54 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _FD_KSTV_H_
#define _FD_KSTV_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

class FDKStV: public FDHookeanMatT, public IsotropicT
{
public:

	/* constructor */
	FDKStV(ifstreamT& in, const FDMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _FD_KSTV_H_ */
