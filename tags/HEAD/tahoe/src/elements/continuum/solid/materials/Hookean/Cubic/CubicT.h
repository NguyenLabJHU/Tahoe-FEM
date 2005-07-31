/* $Id: CubicT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _CUBIC_T_H_
#define _CUBIC_T_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class ifstreamT;
class dMatrixT;

class CubicT
{
public:

	/* constructor */
	CubicT(ifstreamT& in, dMatrixT& moduli);
		
	/* print parameters */
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;

protected:

	/* set modulus */
	void ComputeModuli(dMatrixT& moduli, double C11, double C12, double C44);

protected:

	double fC11;
	double fC12;
	double fC44;
};

#endif /* _CUBIC_T_H_ */
