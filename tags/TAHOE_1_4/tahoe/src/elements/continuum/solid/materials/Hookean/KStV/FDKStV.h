/* $Id: FDKStV.h,v 1.5 2003-01-29 07:34:42 paklein Exp $ */
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
	FDKStV(ifstreamT& in, const FSMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _FD_KSTV_H_ */
