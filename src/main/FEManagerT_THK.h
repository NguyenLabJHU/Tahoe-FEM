/* $Id: FEManagerT_THK.h,v 1.4 2003-05-21 23:47:57 paklein Exp $ */
#ifndef _FE_MANAGER_THK_H_
#define _FE_MANAGER_THK_H_

/* element configuration header */
#include "ElementsConfig.h"
#ifdef BRIDGING_ELEMENT

/* base class */
#include "FEManagerT_bridging.h"
#include "dMatrixT.h"

namespace Tahoe {

/** FEManagerT to support the time history kernel (THK) formulation */
class FEManagerT_THK: public FEManagerT_bridging
{
public:

	/** constructor */
	FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
		ifstreamT& bridging_input);

	/** initialize members */
	virtual void Initialize(InitCodeT init = kFull);

	/** \name solution steps. 
	 * See FEManagerT for more information */
	/*@{*/
	/** initialize the current time increment for all groups */
	virtual ExceptionT::CodeT InitStep(void);

	/** close the current time increment for all groups */
	virtual ExceptionT::CodeT CloseStep(void);
	/*@}*/

	/** contains hard coded values of theta matrix **/
	const dMatrixT& GetTheta(int index);

	/** calculate displacement BC for given THK boundary atoms **/
	const dArray2DT& ComputeStaticDispBC(void);
                                           
private:

        dMatrixT fTheta;
        dArray2DT fBdisp;
        int fNcrit, ncrit;
        iArray2DT fNeighbor0, fNeighbor1;
        iArrayT fNodes, fN0;

};

} /* namespace Tahoe */

#endif /* BRIDGING_ELEMENT */
#endif /* _FE_MANAGER_THK_H_ */
