/* $Id: FEManagerT_THK.cpp,v 1.5 2003-05-05 01:13:44 paklein Exp $ */
#include "FEManagerT_THK.h"
#ifdef BRIDGING_ELEMENT

using namespace Tahoe;

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input)
{

}

#endif /* BRIDGING_ELEMENT */
