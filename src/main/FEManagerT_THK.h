/* $Id: FEManagerT_THK.h,v 1.2 2003-04-05 18:54:37 hspark Exp $ */
#ifndef _FE_MANAGER_THK_H_
#define _FE_MANAGER_THK_H_

/* element configuration header */
#include "ElementsConfig.h"
#ifdef BRIDGING_ELEMENT

/* base class */
#include "FEManagerT_bridging.h"

namespace Tahoe {

/** FEManagerT to support the time history kernel (THK) formulation */
class FEManagerT_THK: public FEManagerT_bridging
{
public:

	/** constructor */
	FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
		ifstreamT& bridging_input);


	/** contains hard coded values of theta matrix **/
	const dMatrixT& GetTheta(int index);

	/** calculate displacement BC for given boundary atom **/
	const dArrayT& ComputeStaticDispBC(const dArray2DT& disp0, const dArray2DT& disp1,
					   int ncrit);

};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_THK_H_ */
