/* $Id: CubicT.h,v 1.4.40.1 2004-01-21 19:10:06 paklein Exp $ */
/* created: paklein (06/11/1997) */
#ifndef _CUBIC_T_H_
#define _CUBIC_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "Material2DT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dMatrixT;

class CubicT: virtual public ParameterInterfaceT
{
public:

	/** constructor */
	CubicT(ifstreamT& in);
	CubicT(void);
		
	/* print parameters */
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

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
