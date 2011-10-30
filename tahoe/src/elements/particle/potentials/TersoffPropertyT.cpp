/* $Id: TersoffPropertyT.cpp,v 1.2.4.1 2011-10-30 06:26:11 bcyansfn Exp $ */
#include "TersoffPropertyT.h"
#include <cstddef>

using namespace Tahoe;

namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<TersoffPropertyT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<TersoffPropertyT>::fByteCopy = false; 
}

/* constructor */
TersoffPropertyT::TersoffPropertyT(void)
{
	SetName("tersoff_property");
}
