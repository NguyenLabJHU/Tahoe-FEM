/* $Id: KBC_PrescribedT.cpp,v 1.3 2004-07-15 08:31:21 paklein Exp $ */
#include "KBC_PrescribedT.h"

using namespace Tahoe;

/* constructor */
KBC_PrescribedT::KBC_PrescribedT(const BasicSupportT& support):
	KBC_ControllerT(support)
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
