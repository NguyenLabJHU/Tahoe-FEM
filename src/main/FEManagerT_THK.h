/* $Id: FEManagerT_THK.h,v 1.3 2003-04-07 06:13:05 paklein Exp $ */
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
};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_THK_H_ */
