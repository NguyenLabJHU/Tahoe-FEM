/* $Id: SSHookean1D.h,v 1.6 2004-07-15 08:27:01 paklein Exp $ */
#ifndef _SS_HOOKEAN_1D_H_
#define _SS_HOOKEAN_1D_H_

/* base classes */
#include "IsotropicT.h"
#include "SSHookeanMatT.h"

namespace Tahoe {

class SSHookean1D: public SSHookeanMatT, public IsotropicT
{
public:

	/** constructor */
	SSHookean1D(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

}//namespace Tahoe

#endif /* _SS_HOOKEAN_1D_H_ */
