/* $Id: ios_fwd_decl.h,v 1.9 2003-11-10 22:14:13 cjkimme Exp $ */
/* created: paklein (08/11/1999) */
/* Include this header instead of writing forward declarations */
/* explicitly. Some compilers do not allow forward declarations */
/* of stream classes. */
#ifndef _IOSFWD_H_
#define _IOSFWD_H_

#include "Environment.h"

#if defined(__JANUS__) && defined(__ROGUE_STL__)
#include <iosfwd>
using namespace std; // can put this here now that namespace std is added
#else
#ifdef _MW_MSL_ // Metrowerks Standard Library
#include <iosfwd.h> //MSL C++ header
#elif defined(__SUNPRO_CC) 
// SUNWspro 5.0
#include <iostream.h>
#include <fstream.h>
#elif defined(__GNU__) && defined (__PGI__)
//PAK(02/28/2001): this header does not seem to work
//#include <iosfwd.h>
#include <iostream.h>
#elif defined(__DEC__) && defined (__USE_STD_IOSTREAM)
#include <iosfwd>
#elif defined(__GCC_3__)
#include <iosfwd>
#elif defined(__INTEL_CC__)
#include <iosfwd>
#else // plain forward declarations OK
class istream;
class ostream;
class ifstream;
class ofstream;
#endif // _MW_MSL_
#endif // __JANUS__ && __ROGUE_STL__
#endif // _IOSFWD_H_
