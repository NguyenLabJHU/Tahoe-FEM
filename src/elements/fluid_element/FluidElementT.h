/* $Header: /home/regueiro/tahoe_cloudforge_repo_snapshots/development/src/elements/fluid_element/FluidElementT.h,v 1.1 2006-07-06 16:28:45 a-kopacz Exp $ */

#ifndef _FLUID_ELEMENT_H_
#define _FLUID_ELEMENT_H_

/* base class */
#include "ContinuumElementT.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class StringT;

/** fluid element */
class FluidElementT: public ContinuumElementT
{
public:

	/* constructor */
	FluidElementT(const ElementSupportT& support);

	/* destructor */
	virtual ~FluidElementT(void);

protected: /* for derived classes only */

	/* called by FormRHS and FormLHS */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);
	virtual void RHSDriver(void);

private:


};

} // namespace Tahoe
#endif /* _FLUID_ELEMENT_H_ */
