/* $Id: HookeanMatT.h,v 1.2 2001-07-03 01:35:05 paklein Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Base class for all Hookean materials, defined as:                      */
/* 	stress_ij = moduli_ijkl strain_kl                                     */

#ifndef _HOOKEAN_MAT_H_
#define _HOOKEAN_MAT_H_

/* forward declarations */
class dSymMatrixT;

/* direct members */
#include "dMatrixT.h"

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

#endif /* _HOOKEAN_MAT_H_ */
