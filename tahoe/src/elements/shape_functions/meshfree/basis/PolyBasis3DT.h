/* $Id: PolyBasis3DT.h,v 1.4 2004-11-03 01:21:01 raregue Exp $ */
/* created: paklein (04/19/2000)                                          */

#ifndef _POLYBASIS_3D_T_H_
#define _POLYBASIS_3D_T_H_

/* base class */
#include "BasisT.h"


namespace Tahoe {

class PolyBasis3DT: public BasisT
{
public:

	/* constructor */
	PolyBasis3DT(int complete);
	
	/* return the number of basis functions */
	virtual int BasisDimension(void) const;

	/* evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

};

} // namespace Tahoe 
#endif /* _POLYBASIS_3D_T_H_ */
