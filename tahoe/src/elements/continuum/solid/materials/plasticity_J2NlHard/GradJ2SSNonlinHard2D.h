/* $Id: GradJ2SSNonlinHard2D.h,v 1.4 2004-09-10 22:39:38 paklein Exp $ */
#ifndef _GRAD_J2_SS_NONLIN_HARD_2D_H_
#define _GRAD_J2_SS_NONLIN_HARD_2D_H_

/* base classes */
#include "GradJ2SSNonlinHard.h"

namespace Tahoe {

class GradJ2SSNonlinHard2D : public GradJ2SSNonlinHard
{
public:

	/* constructor */
	GradJ2SSNonlinHard2D(ifstreamT& in, const SSMatSupportT& support);

	/* returns elastic strain (3D) */
	virtual const dSymMatrixT& ElasticStrain(const dSymMatrixT& totalstrain,
		const ElementCardT& element, int ip);
	
	/* modulus */
	virtual const dMatrixT& c_ijkl(void);
	
	/* stress */
	virtual const dSymMatrixT& s_ij(void);

private:

	/* return values */
	dSymMatrixT fStress2D;
	dMatrixT    fModulus2D;

	/* work space */
	dSymMatrixT	fTotalStrain3D;
};

} // namespace Tahoe 
#endif /* _GRAD_J2_SS_NONLIN_HARD_2D_H_ */
