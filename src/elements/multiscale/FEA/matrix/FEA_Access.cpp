/* $Id: FEA_Access.cpp,v 1.11 2004-11-12 21:00:22 paklein Exp $ */
/** This file contains global parameters for the FEA classes */
#include "FEA.h"

using namespace Tahoe;

#if defined (__DEC__) || (defined (__SUN__) && !defined(__GNUC__)) || defined(__MWERKS__) || defined(__AIX__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT* fStack = NULL;
}

#else

FEA_StackT* fStack = NULL;

#endif
