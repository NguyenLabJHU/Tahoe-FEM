/* $Id: FEManagerT_THK.cpp,v 1.2 2003-04-05 18:54:37 hspark Exp $ */
#include "FEManagerT_THK.h"
#ifdef BRIDGING_ELEMENT

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input)
{

}

/* Obtain theta */
const dMatrixT& FEManagerT_THK::GetTheta(int index)
{





}

/* calculate boundary displacement utilizing theta */
const dArrayT& FEManagerT_THK::ComputeStaticDispBC(const dArray2DT& disp0, const dArray2DT& disp1,
						   int ncrit)
{





}

#endif BRIDGING_ELEMENT
