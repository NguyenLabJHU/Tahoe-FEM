/* $Id: SimoIso3D.h,v 1.10 2004-07-15 08:27:35 paklein Exp $ */
/* created: paklein (03/02/1997) */
#ifndef _SIMO_ISO_3D_H_
#define _SIMO_ISO_3D_H_

/* base classes */
#include "FSIsotropicMatT.h"

namespace Tahoe {

/** hyperelastic material governed by Simo's split volumetric/deviatoric
 * stored energy function.
 * \note This material is inherently 3D
 */
class SimoIso3D: public FSIsotropicMatT
{
public:

	/* constructor */
	SimoIso3D(void);
	
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
//TEMP - no reason to use these in total Lagrangian formulation.
//       calls to these write error message and throw ExceptionT::xception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:
	
	/* computation routines - split volumetric/deviatoric */
	void ComputeModuli(double J, const dSymMatrixT& b_bar, dMatrixT& moduli);
	void ComputeCauchy(double J, const dSymMatrixT& b_bar, dSymMatrixT& cauchy);
	double ComputeEnergy(double J, const dSymMatrixT& b_bar);

	/* Volumetric energy function and derivatives */
	double   U(double J) const;
	double  dU(double J) const;
	double ddU(double J) const;

private:

	/** return true if material implementation supports imposed thermal
	 * strains. This material does support multiplicative thermal
	 * strains. */
	virtual bool SupportsThermalStrain(void) const { return true; };

protected:

	/* return values */
	dSymMatrixT	fStress;
	dMatrixT    fModulus;

	/* work space */
	dSymMatrixT	fb;	
	dSymMatrixT	fb_bar;

private:

	dMatrixT	frank4;
	
	/* fixed forms */
	dSymMatrixT	fIdentity;
	dMatrixT	fIcrossI;
	dMatrixT	fIdentity4;
	dMatrixT	fDevOp4;
};

/* inlines */
inline double SimoIso3D::U(double J) const
{
	return 0.5*Kappa()*(0.5*(J*J - 1.0) - log(J));
}

inline double SimoIso3D::dU(double J) const
{
	return 0.5*Kappa()*(J - 1.0/J);
}

inline double SimoIso3D::ddU(double J) const
{
	return 0.5*Kappa()*(1.0 + 1.0/(J*J));
}

} // namespace Tahoe 
#endif /* _SIMO_ISO_3D_H_ */
