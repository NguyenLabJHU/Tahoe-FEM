/* $Id: ArraySettings.cpp,v 1.9 2003-11-03 18:51:13 paklein Exp $ */
/* created: paklein (01/23/2001) */
#include "ArrayT.h"

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

} /* namespace Tahoe */

#include "RaggedArray2DT.h"

namespace Tahoe {

template<> const bool ArrayT<const RaggedArray2DT<int>*>::fByteCopy = true;

} /* namespace Tahoe */

