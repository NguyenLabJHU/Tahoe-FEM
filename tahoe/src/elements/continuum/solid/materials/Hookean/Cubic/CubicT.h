/* $Id: CubicT.h,v 1.4.40.2 2004-03-02 17:46:13 paklein Exp $ */
/* created: paklein (06/11/1997) */
#ifndef _CUBIC_T_H_
#define _CUBIC_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "SolidMaterialT.h"

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
	void ComputeModuli2D(dMatrixT& moduli, SolidMaterialT::ConstraintT constraint) const;

	/* scale factor for constrained dilatation */
	double DilatationFactor2D(SolidMaterialT::ConstraintT constraint) const;   	

protected:

	double fC11;
	double fC12;
	double fC44;
};

} // namespace Tahoe 
#endif /* _CUBIC_T_H_ */
