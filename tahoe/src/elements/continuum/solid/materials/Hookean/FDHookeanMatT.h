/* $Id: FDHookeanMatT.h,v 1.7.46.1 2004-04-08 07:32:43 paklein Exp $ */
/* created: paklein (06/10/1997) */
#ifndef _FD_HOOKEAN_MAT_H_
#define _FD_HOOKEAN_MAT_H_

/* base classes */
#include "FSSolidMatT.h"
#include "HookeanMatT.h"

namespace Tahoe {

class FDHookeanMatT: public FSSolidMatT, public HookeanMatT
{
public:

	/** constructor */
	FDHookeanMatT(ifstreamT& in, const FSMatSupportT& support);
	FDHookeanMatT(void);

	/** initialization */
	virtual void Initialize(void);

	/** set the material support or pass NULL to clear */
	virtual void SetFSMatSupport(const FSMatSupportT* support);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const;
	/*@}*/

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/* returns the strain energy density for the specified strain */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** return true if material implementation supports imposed thermal
	 * strains. This material does support multiplicative thermal
	 * strains. */
	virtual bool SupportsThermalStrain(void) const { return true; };
	
private:

	/** Green-Lagrangian strain */
	dSymMatrixT fE;

	/** \name return values */
	/*@{*/
	dSymMatrixT fStress;	
	dMatrixT    fModulus;
	FrameT      fLastCall;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _FD_HOOKEAN_MAT_H_ */
