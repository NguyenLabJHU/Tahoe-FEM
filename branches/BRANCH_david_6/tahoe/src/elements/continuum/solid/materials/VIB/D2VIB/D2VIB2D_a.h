/* $Id: D2VIB2D_a.h,v 1.5 2005-12-19 17:16:18 tdnguye Exp $ */
/* created: paklein (10/23/1999) */
#ifndef _D2_VIB_2D_A_H_
#define _D2_VIB_2D_A_H_

/* base class */
#include "D2VIB2D.h"

namespace Tahoe {

class D2VIB2D_a: public D2VIB2D
{
public:

	/* constructor */
	D2VIB2D_a(ifstreamT& in, const D2FSMatSupportT& support);

	/* material internal stress terms */
	virtual void StressTerms(dMatrixT& DW, dMatrixT& DDW);

protected:

	/* nodal displacements */
	const LocalArrayT& fLocDisp;

	/* gradient term coefficient (constant for now) */
	double fD2coeff;

	/* work space */
	dSymMatrixT fPK2;
	dMatrixT fPK2mat;
	dMatrixT fPK1;
	dMatrixT fGradGradU;
};

} // namespace Tahoe 
#endif /* _D2_VIB_2D_A_H_ */
