/* $Id: DOFElementT.cpp,v 1.2 2003-01-27 07:00:24 paklein Exp $ */
/* created: paklein (06/01/1998)                                          */
/* base class to defines the interface for augmented Lagrangian           */
/* element classes (to be used by the corresponding NodeManagerT)         */

#include "DOFElementT.h"
#include "ArrayT.h"

/* array behavior */

using namespace Tahoe;

namespace Tahoe {
const bool ArrayT<DOFElementT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
DOFElementT::DOFElementT(void) { }
