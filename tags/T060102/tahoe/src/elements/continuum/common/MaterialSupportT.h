/* $Id: MaterialSupportT.h,v 1.1 2002-04-21 07:18:27 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _SS_HOOKEAN_MAT_H_
#define _SS_HOOKEAN_MAT_H_

/* base classes */
#include "SSStructMatT.h"
#include "HookeanMatT.h"

class SSHookeanMatT: public SSStructMatT, public HookeanMatT
{
public:

	/* constructor */
	SSHookeanMatT(ifstreamT& in, const SmallStrainT& element);

	/* initialization */
	virtual void Initialize(void);

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
		// overloaded to avoid additional virtual function call

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

protected:

	/* return values */
	dSymMatrixT fStress;
};

#endif /* _SS_HOOKEAN_MAT_H_ */
