/* $Id: FDCubicT.h,v 1.5.46.3 2004-06-09 23:17:35 paklein Exp $ */
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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _FD_CUBIC_T_H_ */
