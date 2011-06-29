/* $Id: DOFElementT.cpp,v 1.4 2003-11-04 01:35:12 paklein Exp $ */
/* created: paklein (06/01/1998)                                          */
/* base class to defines the interface for augmented Lagrangian           */
/* element classes (to be used by the corresponding NodeManagerT)         */

#include "DOFElementT.h"
#include "ArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<DOFElementT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
DOFElementT::DOFElementT(void) { }
