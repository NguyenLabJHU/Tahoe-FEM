/* $Id: ArraySettings.cpp,v 1.7 2002-11-09 01:54:40 paklein Exp $ */
/* created: paklein (01/23/2001) */
#include "ArrayT.h"
#include "RaggedArray2DT.h"

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
const bool ArrayT<const RaggedArray2DT<int>*>::fByteCopy = true;

} /* namespace Tahoe */

