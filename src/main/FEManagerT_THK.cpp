/* $Id: FEManagerT_THK.cpp,v 1.4 2003-05-05 00:14:00 paklein Exp $ */
#include "FEManagerT_THK.h"
#ifdef BRIDGING_ELEMENT

using namespace Tahoe;

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input)
{

}

#endif BRIDGING_ELEMENT
