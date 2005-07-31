/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.1.1.1 2001-07-19 06:24:00 paklein Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#if defined(__STDC__) || defined(__MWERKS__)
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#if 0
/* PAK: this are not system headers */
#include <defs.h>
#include <struct.h>
#include <macros.h>
#include <rename.h>
#include <proto.h>
#endif

#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "rename.h"
#include "proto.h"

#ifdef __cplusplus
}
#endif
