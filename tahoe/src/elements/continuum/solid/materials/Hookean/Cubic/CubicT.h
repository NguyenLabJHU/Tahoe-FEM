/* $Id: CubicT.h,v 1.4 2002-07-05 22:28:15 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */

#ifndef _CUBIC_T_H_
#define _CUBIC_T_H_

#include "Environment.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dMatrixT;

}

/* direct members */
#include "Material2DT.h"

namespace Tahoe {

class CubicT
{
public:

	/* constructor */
	CubicT(ifstreamT& in);
		
	/* print parameters */
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;

protected:

	/* set modulus */
	void ComputeModuli(dMatrixT& moduli);
	void ComputeModuli2D(dMatrixT& moduli, Material2DT::ConstraintOptionT constraint) const;

	/* scale factor for constrained dilatation */
	double DilatationFactor2D(Material2DT::ConstraintOptionT constraint) const;   	

protected:

	double fC11;
	double fC12;
	double fC44;
};

} // namespace Tahoe 
#endif /* _CUBIC_T_H_ */
