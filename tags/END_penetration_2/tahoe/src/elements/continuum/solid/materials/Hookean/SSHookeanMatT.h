/* $Id: SSHookeanMatT.h,v 1.7 2004-01-10 04:41:12 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _SS_HOOKEAN_MAT_H_
#define _SS_HOOKEAN_MAT_H_

/* base classes */
#include "SSSolidMatT.h"
#include "HookeanMatT.h"

namespace Tahoe {

class SSHookeanMatT: public SSSolidMatT, public HookeanMatT
{
public:

	/** constructor */
	SSHookeanMatT(ifstreamT& in, const SSMatSupportT& support);
	SSHookeanMatT(void);

	/** initialization */
	virtual void Initialize(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
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
