/* $Id: FEManagerT_THK.cpp,v 1.3 2003-04-07 06:13:05 paklein Exp $ */
#include "FEManagerT_THK.h"
#ifdef BRIDGING_ELEMENT

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input)
{

}

#endif BRIDGING_ELEMENT
