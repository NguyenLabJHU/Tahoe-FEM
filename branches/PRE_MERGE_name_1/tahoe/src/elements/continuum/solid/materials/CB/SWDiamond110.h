/* $Id: SWDiamond110.h,v 1.2 2001-07-03 01:35:00 paklein Exp $ */
/* created: paklein (08/25/1996)                                          */

#ifndef _SWDIAMOND110_H_
#define _SWDIAMOND110_H_

/* base class */
#include "SWMaterial2D.h"

class SWDiamond110: public SWMaterial2D
{
public:

	/* constructor */
	SWDiamond110(ifstreamT& in, const FiniteStrainT& element);

	/* print name */
	virtual void PrintName(ostream& out) const;
	
protected:	
	
	/* smymetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					
	/* strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

};

#endif /* _SWDIAMOND110_H_ */
