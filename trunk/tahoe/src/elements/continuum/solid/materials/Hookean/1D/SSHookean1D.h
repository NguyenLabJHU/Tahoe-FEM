#ifndef _SS_HOOKEAN_1D_H_
#define _SS_HOOKEAN_1D_H_

/* base classes */
#include "IsotropicT.h"
#include "SSHookeanMatT.h"

namespace Tahoe {

class SSHookean1D: public SSHookeanMatT, public IsotropicT
{
public:

	/* constructor */
	SSHookean1D(ifstreamT& in, const SmallStrainT& element);

	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* set modulus */
	virtual void SetModulus(dMatrixT& modulus);
};

}//namespace Tahoe

#endif /* _SS_HOOKEAN_1D_H_ */
