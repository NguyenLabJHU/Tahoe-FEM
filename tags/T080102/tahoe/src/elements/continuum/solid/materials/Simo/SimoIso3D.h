/* $Id: SimoIso3D.h,v 1.5 2002-07-02 19:55:50 cjkimme Exp $ */
/* created: paklein (03/02/1997)                                          */
/* Hyperelastic material governed by Simo's split volumetric/deviatoric   */
/* stored energy function.                                                */
/* Note: This material is inherently 3D.                                  */

#ifndef _SIMO_ISO_3D_H_
#define _SIMO_ISO_3D_H_

/* base classes */
#include "FDStructMatT.h"
#include "IsotropicT.h"


namespace Tahoe {

class SimoIso3D: public FDStructMatT, public IsotropicT
{
public:

	/* constructor */
	SimoIso3D(ifstreamT& in, const FiniteStrainT& element);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

	/* spatial description */
	virtual const dMatrixT& c_ijkl(void); // spatial tangent moduli
	virtual const dSymMatrixT& s_ij(void); // Cauchy stress

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
//TEMP - no reason to use these in total Lagrangian formulation.
//       calls to these write error message and throw exception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);
	
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
