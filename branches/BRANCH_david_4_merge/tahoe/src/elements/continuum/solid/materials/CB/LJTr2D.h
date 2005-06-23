/* $Id: LJTr2D.h,v 1.8 2004-07-15 08:26:42 paklein Exp $ */
/* created: paklein (07/01/1996) */
#ifndef _LJTR2D_H_
#define _LJTR2D_H_

/* base class */
#include "NL_E_MatT.h"
#include "CBLatticeT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** plane stress hexagonal lattice with LJ potential */
  class LJTr2D: public NL_E_MatT, protected CBLatticeT
{
public:

	/** constructor */
	LJTr2D(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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

	/** LJ scaling constant */
	double feps;
				
};

} // namespace Tahoe 
#endif /* _LJTR2D_H_ */
