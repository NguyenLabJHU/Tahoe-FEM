/* $Id: PolyBasis2DT.h,v 1.3 2004-10-30 20:54:28 raregue Exp $ */
/* created: paklein (12/13/1999)                                          */

#ifndef _POLYBASIS_2D_T_H_
#define _POLYBASIS_2D_T_H_

/* base class */
#include "BasisT.h"


namespace Tahoe {

class PolyBasis2DT: public BasisT
{
public:

	/* constructor */
	PolyBasis2DT(int complete);
	
	/* return the number of basis functions */
	virtual int BasisDimension(void) const;

	/* evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

};

} // namespace Tahoe 
#endif /* _POLYBASIS_2D_T_H_ */
