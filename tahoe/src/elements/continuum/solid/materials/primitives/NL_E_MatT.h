/* $Id: NL_E_MatT.h,v 1.1.1.1.2.2 2001-06-22 14:18:30 paklein Exp $ */
/* created: paklein (06/13/1997)                                          */
/* Base class for materials with nonlinear elastic behavior               */
/* which is computed from Langrangian coordinates (by the pure            */
/* virtual functions below).                                              */
/* Note: The material tangent moduli and 2nd PK are transformed           */
/* to the spatial quantities, spatial tangent moduli and                  */
/* Cauchy stress.                                                         */
/* Note: The moduli, stress, and stored energy are assumed to be          */
/* _at_least_ functions of the Langrangian strain, to ensure that         */
/* the finite deformation continuum is current for the                    */
/* required material->spatial transformations. Any other                  */
/* dependencies, ie. visco-elasticity must be handled by                  */
/* the derived material classes.                                          */
/* Note: The particular material orientation with respect to the          */
/* global axes is also assumed to be taken care of by the                 */
/* derived material classes.                                              */

#ifndef _NL_E_MAT_T_H_
#define _NL_E_MAT_T_H_

/* base classes */
#include "FDStructMatT.h"

class NL_E_MatT: public FDStructMatT
{
  public:

	/** constructor */
	NL_E_MatT(ifstreamT& in, const FiniteStrainT& element);
	
	/* spatial description */
	virtual const dMatrixT& c_ijkl(void);  /**< spatial tangent modulus */
	virtual const dSymMatrixT& s_ij(void); /**< Cauchy stress */

	/* material description */
	virtual const dMatrixT& C_IJKL(void);  /**< material tangent moduli */
	virtual const dSymMatrixT& S_IJ(void); /**< 2nd Piola-Kirchhoff stress */

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

#endif /* _NL_E_MAT_T_H_ */
