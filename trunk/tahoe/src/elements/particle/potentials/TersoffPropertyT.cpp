/* $Id: TersoffPropertyT.cpp,v 1.3 2011-12-01 21:11:39 bcyansfn Exp $ */
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
