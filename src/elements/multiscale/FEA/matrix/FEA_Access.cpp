/* $Id: FEA_Access.cpp,v 1.9 2003-04-11 23:05:54 paklein Exp $ */
/** This file contains global parameters for the FEA classes */
#include "FEA.h"

#if defined (__DEC__) || (defined (__SUN__) && !defined(__GNU__)) || defined(__MWERKS__) || defined(__AIX__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT* fStack = NULL;
}

#else

FEA_StackT* fStack = NULL;

#endif
