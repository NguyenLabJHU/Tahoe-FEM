/* $Id: DOFElementT.cpp,v 1.3 2003-10-28 07:28:09 paklein Exp $ */
/* created: paklein (06/01/1998)                                          */
/* base class to defines the interface for augmented Lagrangian           */
/* element classes (to be used by the corresponding NodeManagerT)         */

#include "DOFElementT.h"
#include "ArrayT.h"

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
template<> const bool ArrayT<DOFElementT*>::fByteCopy = true;
} /* namespace Tahoe */

/* constructor */
DOFElementT::DOFElementT(void) { }
