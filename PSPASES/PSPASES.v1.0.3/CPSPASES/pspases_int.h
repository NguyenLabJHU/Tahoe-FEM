/* $Id: pspases_int.h,v 1.1 2005-01-05 07:37:08 paklein Exp $ */

#ifndef PSPASES_INT_H
#define PSPASES_INT_H

#include "pspases_f2c.h"

/* not debugging */
#ifdef NDEBUG

/* skip wrappers */
#define mydsyrk_ dsyrk_

#endif

#endif /* PSPASES_INT_H */
