/* $Id: SWDiamond100.h,v 1.3 2002-07-02 19:55:34 cjkimme Exp $ */
/* created: paklein (08/25/1996)                                          */

#ifndef _SWDIAMOND100_H_
#define _SWDIAMOND100_H_

/* base class */
#include "SWMaterial2D.h"


namespace Tahoe {

class SWDiamond100: public SWMaterial2D
{
public:

	/* constructor */
	SWDiamond100(ifstreamT& in, const FiniteStrainT& element);

	/* print name */
	virtual void PrintName(ostream& out) const;

protected:	
	
	/* symmetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					
	/* strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

};

} // namespace Tahoe 
#endif /* _SWDIAMOND100_H_ */
