/* $Id: FEManagerT_THK.cpp,v 1.1 2003-04-05 08:38:38 paklein Exp $ */
#include "FEManagerT_THK.h"
#ifdef BRIDGING_ELEMENT

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input)
{

}

#endif BRIDGING_ELEMENT
