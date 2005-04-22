/* $Id: PolyBasis2DT.h,v 1.6 2005-04-22 00:33:40 paklein Exp $ */
/* created: paklein (12/13/1999) */
#ifndef _POLYBASIS_2D_T_H_
#define _POLYBASIS_2D_T_H_

/* base class */
#include "BasisT.h"

namespace Tahoe {

/** monomial basis in two dimensions */
class PolyBasis2DT: public BasisT
{
public:

	/** constructor */
	PolyBasis2DT(int complete, bool cross_terms);
	
	/** return the number of basis functions */
	virtual int BasisDimension(void) const;

	/** evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

private:

	/** include monomials for cross terms */
	bool fCrossTerms;
};

} // namespace Tahoe 
#endif /* _POLYBASIS_2D_T_H_ */
