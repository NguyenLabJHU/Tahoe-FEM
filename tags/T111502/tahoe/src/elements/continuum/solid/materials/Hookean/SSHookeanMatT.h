/* $Id: SSHookeanMatT.h,v 1.5 2002-11-14 17:06:02 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _SS_HOOKEAN_MAT_H_
#define _SS_HOOKEAN_MAT_H_

/* base classes */
#include "SSStructMatT.h"
#include "HookeanMatT.h"

namespace Tahoe {

class SSHookeanMatT: public SSStructMatT, public HookeanMatT
{
public:

	/* constructor */
	SSHookeanMatT(ifstreamT& in, const SSMatSupportT& support);

	/* initialization */
	virtual void Initialize(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * StructuralMaterialT::s_ij. See StructuralMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fStress.Trace()/3.0; };
	/*@}*/

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

} // namespace Tahoe 
#endif /* _SS_HOOKEAN_MAT_H_ */
