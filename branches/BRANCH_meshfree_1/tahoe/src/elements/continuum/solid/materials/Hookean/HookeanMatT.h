/* $Id: HookeanMatT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Base class for all Hookean materials, defined as:                      */
/* 	stress_ij = moduli_ijkl strain_kl                                     */

#ifndef _HOOKEAN_MAT_H_
#define _HOOKEAN_MAT_H_

/* forward declarations */
class dMatrixT;
class dSymMatrixT;

class HookeanMatT
{
public:

	/* constructor */
	HookeanMatT(void);

protected:

	/* symmetric stress */
	void HookeanStress(const dMatrixT& moduli, const dSymMatrixT& strain,
		dSymMatrixT& stress);

	/* strain energy density */
	double HookeanEnergy(const dMatrixT& moduli,
		const dSymMatrixT& strain);
};

#endif /* _HOOKEAN_MAT_H_ */
