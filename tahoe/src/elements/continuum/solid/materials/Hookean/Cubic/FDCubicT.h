/* $Id: FDCubicT.h,v 1.5.30.1 2004-01-21 19:10:06 paklein Exp $ */
/* created: paklein (06/11/1997) */
#ifndef _FD_CUBIC_T_H_
#define _FD_CUBIC_T_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "CubicT.h"

namespace Tahoe {

class FDCubicT: public FDHookeanMatT, public CubicT
{
public:

	/** constructor */
	FDCubicT(ifstreamT& in, const FSMatSupportT& support);
	FDCubicT(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _FD_CUBIC_T_H_ */
