/* $Id: SSKStV.h,v 1.2 2001-07-03 01:35:09 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _SS_KSTV_H_
#define _SS_KSTV_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "IsotropicT.h"

class SSKStV: public SSHookeanMatT, public IsotropicT
{
public:

	/* constructor */
	SSKStV(ifstreamT& in, const SmallStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

#endif /* _SS_KSTV_H_ */
