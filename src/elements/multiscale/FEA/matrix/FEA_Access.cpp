/* $Id: FEA_Access.cpp,v 1.6 2002-12-02 07:11:03 paklein Exp $ */
/** This file contains global parameters for the FEA classes */
//DEVELOPMENT
#include "FEA.h"

#if defined (__DEC__) || defined (__SUN__) || defined(__MWERKS__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT* fStack = NULL;
}

#else

FEA_StackT* fStack = NULL;

#endif
