// ios_fwd_decl.h
//
// include this header instead of writing forward declarations
// explicitly. MSL does not allow forward declarations of stream
// classes.

// created:       PAK (08/11/1999)
// last modified: PAK (08/11/1999)

#ifndef _IOSFWD_H_
#define _IOSFWD_H_

#include "Environment.h"

#ifdef _MW_MSL_ // Metrowerks Standard Library
#include <iosfwd> //MSL C++ header
class ifstream;
class ofstream;
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
