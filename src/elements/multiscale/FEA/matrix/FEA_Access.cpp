/* $Id: FEA_Access.cpp,v 1.8 2003-02-03 04:40:24 paklein Exp $ */
/** This file contains global parameters for the FEA classes */
#include "FEA.h"

#if defined (__DEC__) || defined (__SUN__) || defined(__MWERKS__) || defined(__AIX__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT* fStack = NULL;
}

#else

FEA_StackT* fStack = NULL;

#endif
