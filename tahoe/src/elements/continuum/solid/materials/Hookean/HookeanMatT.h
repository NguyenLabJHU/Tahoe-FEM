/* $Id: HookeanMatT.h,v 1.5 2004-01-10 04:41:12 paklein Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _HOOKEAN_MAT_H_
#define _HOOKEAN_MAT_H_

/* direct members */
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class dSymMatrixT;

/** base class for all Hookean materials, defined as:
 * 	stress_ij = moduli_ijkl strain_kl */
class HookeanMatT
{
public:

	/** constructor */
	HookeanMatT(int nsd);
	HookeanMatT(void) { };

	/** destructor */
	virtual ~HookeanMatT(void);

	/** initialization */
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
