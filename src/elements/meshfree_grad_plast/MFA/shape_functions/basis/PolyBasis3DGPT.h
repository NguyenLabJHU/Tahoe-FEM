/* $Id: PolyBasis3DGPT.h,v 1.1 2004-06-22 23:17:48 kyonten Exp $ */
/* created: paklein (04/19/2000)                                          */

#ifndef _POLYBASIS_3D_GP_T_H_
#define _POLYBASIS_3D_GP_T_H_

/* base class */
#include "BasisGPT.h"


namespace Tahoe {

class PolyBasis3DGPT: public BasisGPT
{
public:

	/* constructor */
	PolyBasis3DGPT(int complete);
	
	/* return the number of basis functions */
	virtual int BasisDimension(void) const;

	/* evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

};

} // namespace Tahoe 
#endif /* _POLYBASIS_3D_GP_T_H_ */
