/* $Id: SSCubicT.h,v 1.4.48.1 2004-04-08 07:32:46 paklein Exp $ */
/* created: paklein (06/11/1997) */
#ifndef _SS_CUBIC_T_H_
#define _SS_CUBIC_T_H_

/* base classes */
#include "SSHookeanMatT.h"
#include "CubicT.h"

namespace Tahoe {

/** elastic small strain material with cubic symmetry */
class SSCubicT: public SSHookeanMatT, public CubicT
{
public:

	/** constructor */
	SSCubicT(ifstreamT& in, const SSMatSupportT& support);
	SSCubicT(void);

	/** print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _SS_CUBIC_T_H_ */
