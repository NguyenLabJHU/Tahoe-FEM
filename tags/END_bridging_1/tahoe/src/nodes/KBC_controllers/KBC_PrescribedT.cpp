/* $Id: KBC_PrescribedT.cpp,v 1.1.2.1 2003-02-10 02:12:23 paklein Exp $ */
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
