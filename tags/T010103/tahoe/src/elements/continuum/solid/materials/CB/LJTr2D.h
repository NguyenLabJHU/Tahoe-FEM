/* $Id: LJTr2D.h,v 1.6 2002-12-18 22:49:59 cjkimme Exp $ */
/* created: paklein (07/01/1996) */
#ifndef _LJTR2D_H_
#define _LJTR2D_H_

/* base class */
#include "NL_E_Mat2DT.h"
#include "CBLatticeT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** plane stress hexagonal lattice with LJ potential */
  class LJTr2D: public NL_E_Mat2DT, protected CBLatticeT
{
public:

	/* constructor */
	LJTr2D(ifstreamT& in, const FDMatSupportT& support);
	
	virtual void Initialize(void);

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

	virtual void LoadBondTable(void);

private:

	double Ulj(double r) const;
	
	double dUlj(double r) const;

	/* second derivative of the Lennard-Jones 6/12 potential */
	double ddU(double l) const;

private:

	/* LJ scaling constant */
	double feps;
				
};

} // namespace Tahoe 
#endif /* _LJTR2D_H_ */
