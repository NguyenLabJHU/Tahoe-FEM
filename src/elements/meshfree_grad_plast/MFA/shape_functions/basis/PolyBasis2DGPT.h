/* $Id: PolyBasis2DGPT.h,v 1.1 2004-06-22 23:17:48 kyonten Exp $ */
/* created: paklein (12/13/1999)                                          */

#ifndef _POLYBASIS_2D_GP_T_H_
#define _POLYBASIS_2D_GP_T_H_

/* base class */
#include "BasisGPT.h"


namespace Tahoe {

class PolyBasis2DGPT: public BasisGPT
{
public:

	/* constructor */
	PolyBasis2DGPT(int complete);
	
	/* return the number of basis functions */
	virtual int BasisDimension(void) const;

	/* evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

};

} // namespace Tahoe 
#endif /* _POLYBASIS_2D_GP_T_H_ */
