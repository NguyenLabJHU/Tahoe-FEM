/* $Id: FDHookeanMatT.h,v 1.3.6.1 2002-06-27 18:03:09 cjkimme Exp $ */
/* created: paklein (06/10/1997) */

#ifndef _FD_HOOKEAN_MAT_H_
#define _FD_HOOKEAN_MAT_H_

/* base classes */
#include "FDStructMatT.h"
#include "HookeanMatT.h"


namespace Tahoe {

class FDHookeanMatT: public FDStructMatT, public HookeanMatT
{
public:

	/* constructor */
	FDHookeanMatT(ifstreamT& in, const FiniteStrainT& element);

	/* initialization */
	virtual void Initialize(void);

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

private:

	/** return true if material implementation supports imposed thermal
	 * strains. This material does support multiplicative thermal
	 * strains. */
	virtual bool SupportsThermalStrain(void) const { return true; };
	
private:

	/* Green-Lagrangian strain */
	dSymMatrixT fE;

	/* return values */
	dSymMatrixT fStress;	
	dMatrixT    fModulus;	
};

} // namespace Tahoe 
#endif /* _FD_HOOKEAN_MAT_H_ */
