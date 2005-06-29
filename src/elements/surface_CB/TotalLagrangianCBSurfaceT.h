/* $Id: TotalLagrangianCBSurfaceT.h,v 1.1 2005-06-29 17:39:46 paklein Exp $ */
#ifndef _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_
#define _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_

/* base class */
#include "TotalLagrangianT.h"

namespace Tahoe {

/** total Lagrangian, finite strain element for working with Cauchy-Born approach
 * for modeling surface effects */
class TotalLagrangianCBSurfaceT: public TotalLagrangianT
{
public:

	/** constructors */
	TotalLagrangianCBSurfaceT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
		
protected:

	/** list of elements on the surface */
	iArrayT fSurfaceElements;

	/** elements neighbors */
	iArray2DT fSurfaceElementNeighbors;
};

} /* namespace Tahoe */

#endif /* _TOTAL_LAGRANGRIAN_CB_SURFACE_T_H_ */
