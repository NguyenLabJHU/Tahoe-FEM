/* $Id: TotalLagrangianFlatT.h,v 1.2 2004-07-15 08:29:02 paklein Exp $ */
#ifndef _TOTAL_LAGRANGRIAN_FLAT_T_H_
#define _TOTAL_LAGRANGRIAN_FLAT_T_H_

/* base class */
#include "TotalLagrangianT.h"

namespace Tahoe {

/** TotalLagrangianT with a flattened method to compute the residual force */
class TotalLagrangianFlatT: public TotalLagrangianT
{
public:
	
	/** constructor */
	TotalLagrangianFlatT(const ElementSupportT& support);

protected:

	/** form the residual force vector. This is a flattened version of
	 * SolidElementT::RHSDriver combined with SolidElementT::ElementRHSDriver.
	 * The implementation assumes the constitutive routines need information
	 * about the deformation at the current integration point only. */
	virtual void RHSDriver(void);
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_FLAT_T_H_ */
