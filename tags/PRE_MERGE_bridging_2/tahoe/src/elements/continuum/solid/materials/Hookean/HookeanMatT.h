/* $Id: HookeanMatT.h,v 1.4 2002-07-05 22:28:14 paklein Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Base class for all Hookean materials, defined as:                      */
/* 	stress_ij = moduli_ijkl strain_kl                                     */

#ifndef _HOOKEAN_MAT_H_
#define _HOOKEAN_MAT_H_


namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

}

/* direct members */
#include "dMatrixT.h"

namespace Tahoe {

class HookeanMatT
{
public:

	/* constructor */
	HookeanMatT(int nsd);

	/* destructor */
	virtual ~HookeanMatT(void);

	/* initialization */
	void Initialize(void);

protected:

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus) = 0;
	const dMatrixT& Modulus(void) const { return fModulus; };

	/* symmetric stress */
	void HookeanStress(const dSymMatrixT& strain, dSymMatrixT& stress) const;

	/* strain energy density */
	double HookeanEnergy(const dSymMatrixT& strain) const;
		
private:

	/* (constant) modulus */
	dMatrixT fModulus;
};

} // namespace Tahoe 
#endif /* _HOOKEAN_MAT_H_ */
