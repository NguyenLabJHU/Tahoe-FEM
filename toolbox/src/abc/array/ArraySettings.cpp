/* $Id: ArraySettings.cpp,v 1.6 2002-07-05 17:15:48 paklein Exp $ */
/* created: paklein (01/23/2001) */

#include "ArrayT.h"

/* NOTE: IBM's Visual Age C++ compiler xlC requires template
 *       specializations to be declared in the same namespace
 *       as the template. Other platforms don't seem to care. */

namespace Tahoe {
/* built-in types */
const bool ArrayT<int>::fByteCopy = true;
const bool ArrayT<char>::fByteCopy = true;
const bool ArrayT<bool>::fByteCopy = true;
const bool ArrayT<float>::fByteCopy = true;
const bool ArrayT<double>::fByteCopy = true;

/* and their pointers */
const bool ArrayT<int*>::fByteCopy = true;
const bool ArrayT<char*>::fByteCopy = true;
const bool ArrayT<bool*>::fByteCopy = true;
const bool ArrayT<void*>::fByteCopy = true;
const bool ArrayT<float*>::fByteCopy = true;
const bool ArrayT<double*>::fByteCopy = true;

/* arrays of arrays */
const bool ArrayT<ArrayT<int>*>::fByteCopy = true;
const bool ArrayT<ArrayT<double>*>::fByteCopy = true;
} /* namespace Tahoe */

