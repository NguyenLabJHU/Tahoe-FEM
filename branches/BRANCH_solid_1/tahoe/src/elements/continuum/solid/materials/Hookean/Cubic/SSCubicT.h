/* $Id: SSCubicT.h,v 1.1.1.1.2.1 2001-06-06 16:22:02 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _SS_CUBIC_T_H_
#define _SS_CUBIC_T_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "CubicT.h"

class SSCubicT: public SSHookeanMatT, public CubicT
{
public:

	/* constructor */
	SSCubicT(ifstreamT& in, const ElasticT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

#endif /* _SS_CUBIC_T_H_ */
