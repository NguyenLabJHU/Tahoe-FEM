/* $Id: KBC_PrescribedT.cpp,v 1.2 2003-03-31 23:02:50 paklein Exp $ */
#include "KBC_PrescribedT.h"

using namespace Tahoe;

/* constructor */
KBC_PrescribedT::KBC_PrescribedT(NodeManagerT& node_manager):
	KBC_ControllerT(node_manager)
{

}

void KBC_PrescribedT::Initialize(ifstreamT& in)
{
#pragma unused(in)

	// nothing required	
}

void KBC_PrescribedT::InitialCondition(void)
{
	// nothing required	
}
