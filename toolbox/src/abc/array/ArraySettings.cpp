/* $Id: ArraySettings.cpp,v 1.8 2003-10-27 19:50:33 paklein Exp $ */
/* created: paklein (01/23/2001) */
#include "ArrayT.h"
#include "RaggedArray2DT.h"

/* NOTE: IBM's Visual Age C++ compiler xlC requires template
 *       specializations to be declared in the same namespace
 *       as the template. Other platforms don't seem to care. */

namespace Tahoe {

/* built-in types */
template<> const bool ArrayT<int>::fByteCopy = true;
template<> const bool ArrayT<char>::fByteCopy = true;
template<> const bool ArrayT<bool>::fByteCopy = true;
template<> const bool ArrayT<float>::fByteCopy = true;
template<> const bool ArrayT<double>::fByteCopy = true;

/* and their pointers */
template<> const bool ArrayT<int*>::fByteCopy = true;
template<> const bool ArrayT<char*>::fByteCopy = true;
template<> const bool ArrayT<bool*>::fByteCopy = true;
template<> const bool ArrayT<void*>::fByteCopy = true;
template<> const bool ArrayT<float*>::fByteCopy = true;
template<> const bool ArrayT<double*>::fByteCopy = true;

/* arrays of arrays */
template<> const bool ArrayT<ArrayT<int>*>::fByteCopy = true;
template<> const bool ArrayT<ArrayT<double>*>::fByteCopy = true;
template<> const bool ArrayT<const RaggedArray2DT<int>*>::fByteCopy = true;

} /* namespace Tahoe */

