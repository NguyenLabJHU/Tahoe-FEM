/* $Id: FEA_Access.cpp,v 1.7 2002-12-15 01:10:24 paklein Exp $ */
/** This file contains global parameters for the FEA classes */
//DEVELOPMENT
#include "FEA.h"

#if defined (__DEC__) || defined (__SUN__) || defined(__MWERKS__) || defined(__AIX__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT* fStack = NULL;
}

#else

FEA_StackT* fStack = NULL;

#endif
