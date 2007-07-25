/* $Id: TersoffPropertyT.cpp,v 1.2 2006-07-25 16:29:47 d-farrell2 Exp $ */
#include "TersoffPropertyT.h"
#include <stddef.h>

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
