/* $Id: FEA_Access.cpp,v 1.10 2003-05-05 00:58:07 paklein Exp $ */
/** This file contains global parameters for the FEA classes */
#include "FEA.h"

using namespace Tahoe;

#if defined (__DEC__) || (defined (__SUN__) && !defined(__GNU__)) || defined(__MWERKS__) || defined(__AIX__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT* fStack = NULL;
}

#else

FEA_StackT* fStack = NULL;

#endif
