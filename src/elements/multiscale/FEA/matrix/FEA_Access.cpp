/** This file contains global parameters for the FEA classes */

//DEVELOPMENT

#include "FEA.h"

#ifdef __DEC__

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT fStack;
}

#else

FEA_StackT fStack;

#endif
