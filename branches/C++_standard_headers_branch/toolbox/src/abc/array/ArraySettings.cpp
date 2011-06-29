/* $Id: ArraySettings.cpp,v 1.11 2004-04-27 07:23:24 paklein Exp $ */
/* created: paklein (01/23/2001) */
#include "ArrayT.h"

namespace Tahoe {

/* built-in types */
DEFINE_TEMPLATE_STATIC const bool ArrayT<int>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<char>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<bool>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<float>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<double>::fByteCopy = true;

/* and their pointers */
DEFINE_TEMPLATE_STATIC const bool ArrayT<int*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<char*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<bool*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<void*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<float*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<double*>::fByteCopy = true;

/* arrays of arrays */
DEFINE_TEMPLATE_STATIC const bool ArrayT<ArrayT<int>*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<ArrayT<double>*>::fByteCopy = true;

} /* namespace Tahoe */

#include "RaggedArray2DT.h"
#include "nVariArray2DT.h"

namespace Tahoe {

DEFINE_TEMPLATE_STATIC const bool ArrayT<const RaggedArray2DT<int>*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<nVariArray2DT<int>*>::fByteCopy = true;

} /* namespace Tahoe */

