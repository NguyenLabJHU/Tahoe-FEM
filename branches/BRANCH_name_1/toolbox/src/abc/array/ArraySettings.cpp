/* $Id: ArraySettings.cpp,v 1.3.4.1 2002-06-27 18:00:44 cjkimme Exp $ */
/* created: paklein (01/23/2001) */

#include "ArrayT.h"


using namespace Tahoe;

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
