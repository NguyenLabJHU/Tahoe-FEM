/* $Id: NL_E_MatT.h,v 1.6 2003-01-29 07:35:08 paklein Exp $ */
/* created: paklein (06/13/1997) */
#ifndef _NL_E_MAT_T_H_
#define _NL_E_MAT_T_H_

/* base classes */
#include "FSSolidMatT.h"

namespace Tahoe {

/** Base class for materials with nonlinear elastic behavior
 * which is computed from Langrangian coordinates (by the pure
 * virtual functions below).
 * Note: The material tangent moduli and 2nd PK are transformed
 * to the spatial quantities, spatial tangent moduli and
 * Cauchy stress.
 * Note: The moduli, stress, and stored energy are assumed to be
 * _at_least_ functions of the Langrangian strain, to ensure that
 * the finite deformation continuum is current for the
 * required material->spatial transformations. Any other
 * dependencies, ie. visco-elasticity must be handled by
 * the derived material classes.
 * Note: The particular material orientation with respect to the
 * global axes is also assumed to be taken care of by the
 * derived material classes. */
class NL_E_MatT: public FSSolidMatT
{
  public:

	/** constructor */
	NL_E_MatT(ifstreamT& in, const FSMatSupportT& support);
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. \note only correct if NL_E_MatT::s_ij
	 * called previously, otherwise get the trace of S_IJ. */
	virtual double Pressure(void) const { return fPK2.Trace()/3.0; };
	/*@}*/

	/** \name material description */
	/*@{*/
	/** material tangent moduli */
	virtual const dMatrixT& C_IJKL(void);

	/** 2nd Piola-Kirchhoff stress */
	virtual const dSymMatrixT& S_IJ(void);
	/*@}*/

	/** strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli) = 0;	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2) = 0;

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E) = 0;

protected:

	/* Green-Lagrangian strain */
	dSymMatrixT fE;

	/* return values */
	dSymMatrixT	fPK2;
	dMatrixT fModuli;
};

} // namespace Tahoe 
#endif /* _NL_E_MAT_T_H_ */
