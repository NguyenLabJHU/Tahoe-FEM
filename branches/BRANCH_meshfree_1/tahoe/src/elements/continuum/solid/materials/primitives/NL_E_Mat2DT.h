/* $Id: NL_E_Mat2DT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/13/1997)                                          */
/* Base class for materials with 2D nonlinear elastic behavior.           */
/* (See notes in NL_E_MatT.h)                                             */
/* Note: The 2D deformation constraint is assumed to be prescribed        */
/* by the derived material classes, ie. it is not specified               */
/* by an input file parameter.                                            */

#ifndef _NL_E_MAT_2D_T_H_
#define _NL_E_MAT_2D_T_H_

/* base classes */
#include "NL_E_MatT.h"
#include "Material2DT.h"

class NL_E_Mat2DT: public NL_E_MatT, public Material2DT
{
public:

	/* constructor */
	NL_E_Mat2DT(ifstreamT& in, const ElasticT& element, ConstraintOptionT constraint);

	/* print parameters */
	virtual void Print(ostream& out) const;
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

protected:

	/* symmetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli) = 0;	
	
	/* symmetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2) = 0;

	/* strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E) = 0;
	
};

#endif /* _NL_E_MAT_2D_T_H_ */
