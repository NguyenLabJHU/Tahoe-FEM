/* $Id: FDKStV.h,v 1.5.30.2 2004-03-04 06:45:25 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _FD_KSTV_H_
#define _FD_KSTV_H_

/* base classes */
#include "FDHookeanMatT.h"
#include "IsotropicT.h"

namespace Tahoe {

class FDKStV: public FDHookeanMatT, public IsotropicT
{
public:

	/** constructor */
	FDKStV(ifstreamT& in, const FSMatSupportT& support);
	FDKStV(void);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

} // namespace Tahoe 
#endif /* _FD_KSTV_H_ */
