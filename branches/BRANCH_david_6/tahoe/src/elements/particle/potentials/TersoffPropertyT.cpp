/* $Id: TersoffPropertyT.cpp,v 1.1.2.1 2006-06-12 18:40:06 d-farrell2 Exp $ */
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
