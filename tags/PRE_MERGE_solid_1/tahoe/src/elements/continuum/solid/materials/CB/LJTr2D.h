/* $Id: LJTr2D.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (07/01/1996)                                          */
/* Plane stress hexagonal lattice with LJ potential                       */

#ifndef _LJTR2D_H_
#define _LJTR2D_H_

/* base class */
#include "NL_E_RotMat2DT.h"

/* direct members */
#include "dArray2DT.h"

class LJTr2D: public NL_E_RotMat2DT
{
public:

	/* constructor */
	LJTr2D(ifstreamT& in, const ElasticT& element);
	
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

#endif /* _LJTR2D_H_ */
