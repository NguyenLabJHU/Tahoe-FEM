/* $Id: LJTr2D.h,v 1.4 2002-11-14 17:06:00 paklein Exp $ */
/* created: paklein (07/01/1996) */
#ifndef _LJTR2D_H_
#define _LJTR2D_H_

/* base class */
#include "NL_E_RotMat2DT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** plane stress hexagonal lattice with LJ potential */
class LJTr2D: public NL_E_RotMat2DT
{
public:

	/* constructor */
	LJTr2D(ifstreamT& in, const FDMatSupportT& support);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					                 					
	/* strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

private:

	/* second derivative of the Lennard-Jones 6/12 potential */
	double ddU(double l) const;

private:

	/* LJ scaling constant */
	double fScale;
	
	/* bond vectors (undeformed) */
	dArray2DT fBondVectors;
				
};

} // namespace Tahoe 
#endif /* _LJTR2D_H_ */
