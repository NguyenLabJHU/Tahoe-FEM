/** This file contains global parameters for the FEA classes */

//DEVELOPMENT

#include "FEA.h"

#if defined (__DEC__) || defined (__SUN__)

/* declare "global" within the Tahoe namespace */
namespace Tahoe {
FEA_StackT fStack;
}

#else

FEA_StackT fStack;

#endif
