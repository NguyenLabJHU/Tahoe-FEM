/* $Id: SWDiamond110.h,v 1.7.2.1 2004-12-26 06:27:46 d-farrell2 Exp $ */
/* created: paklein (08/25/1996) */
#ifndef _SWDIAMOND110_H_
#define _SWDIAMOND110_H_

/* base class */
#include "SWMaterial2D.h"

namespace Tahoe {

class SWDiamond110: public SWMaterial2D
{
public:

	/* constructor */
	SWDiamond110(ifstreamT& in, const FSMatSupportT& support);
	
protected:	
	
	/* smymetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* symmetric 2nd Piola-Kirchhoff stress */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);
					
	/* strain energy density */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);

};

} // namespace Tahoe 
#endif /* _SWDIAMOND110_H_ */
