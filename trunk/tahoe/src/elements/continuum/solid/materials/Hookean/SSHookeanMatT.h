/* $Id: SSHookeanMatT.h,v 1.8 2004-07-15 08:26:56 paklein Exp $ */
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
	SSHookeanMatT(void);

	/** set the material support or pass NULL to clear */
	virtual void SetSSMatSupport(const SSMatSupportT* support);

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

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* return values */
	dSymMatrixT fStress;
};

} // namespace Tahoe 
#endif /* _SS_HOOKEAN_MAT_H_ */
