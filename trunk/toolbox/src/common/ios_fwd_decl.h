/* $Id: ios_fwd_decl.h,v 1.1.1.1 2001-01-25 20:56:28 paklein Exp $ */
/* created: paklein (08/11/1999)                                          */
/* iosfwd.h                                                               */
/* include this header instead of writing forward declarations            */
/* explicitly. MSL does not allow forward declarations of stream          */
/* classes.                                                               */

#ifndef _IOSFWD_H_
#define _IOSFWD_H_

#include "Environment.h"

#ifdef _MW_MSL_ // Metrowerks Standard Library
#include <iosfwd.h> //MSL C++ header
#else
#ifdef __SUNPRO_CC // SUNWspro 5.0
#include <iostream.h>
#include <fstream.h>
#else // forward declarations OK
class istream;
class ostream;
class ifstream;
class ofstream;
#endif // __SUNPRO_CC
#endif // _MW_MSL_

#endif // _IOSFWD_H_
