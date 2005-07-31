/* $Id: PolyBasis1DT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (12/11/1999)                                          */
/* base class for basis functions                                         */

#ifndef _POLYBASIS_1D_T_H_
#define _POLYBASIS_1D_T_H_

/* base class */
#include "BasisT.h"

class PolyBasis1DT: public BasisT
{
public:

	/* constructor */
	PolyBasis1DT(int complete);
	
	/* return the number of basis functions */
	virtual int BasisDimension(void) const;

	/* evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

};

#endif /* _POLYBASIS_1D_T_H_ */
