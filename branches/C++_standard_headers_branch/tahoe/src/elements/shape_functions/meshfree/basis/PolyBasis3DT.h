/* $Id: PolyBasis3DT.h,v 1.6 2005-04-22 00:33:40 paklein Exp $ */
/* created: paklein (04/19/2000) */
#ifndef _POLYBASIS_3D_T_H_
#define _POLYBASIS_3D_T_H_

/* base class */
#include "BasisT.h"

namespace Tahoe {

/** monomial basis in three dimensions */ 
class PolyBasis3DT: public BasisT
{
public:

	/** constructor */
	PolyBasis3DT(int complete, bool cross_terms);
	
	/** return the number of basis functions */
	virtual int BasisDimension(void) const;

	/** evaluate basis functions at coords */
	virtual void SetBasis(const dArray2DT& coords, int order);

private:

	/** include monomials for cross terms */
	bool fCrossTerms;

};

} // namespace Tahoe 
#endif /* _POLYBASIS_3D_T_H_ */
