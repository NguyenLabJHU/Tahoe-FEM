/* $Id: SSKStV.h,v 1.3.8.1 2002-10-28 06:48:54 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _SS_KSTV_H_
#define _SS_KSTV_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

class SSKStV: public SSHookeanMatT, public IsotropicT
{
public:

	/* constructor */
	SSKStV(ifstreamT& in, const SSMatSupportT& support);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _SS_KSTV_H_ */
