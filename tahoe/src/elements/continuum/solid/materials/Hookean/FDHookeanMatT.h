/* $Id: FDHookeanMatT.h,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/10/1997)                                          */

#ifndef _FD_HOOKEAN_MAT_H_
#define _FD_HOOKEAN_MAT_H_

/* base classes */
#include "FDStructMatT.h"
#include "HookeanMatT.h"

class FDHookeanMatT: public FDStructMatT, public HookeanMatT
{
public:

	/* constructor */
	FDHookeanMatT(ifstreamT& in, const ElasticT& element);

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);
	
protected:

	/* return values */
	dSymMatrixT fStress;
	dMatrixT    fModulus;	
};

#endif /* _FD_HOOKEAN_MAT_H_ */
