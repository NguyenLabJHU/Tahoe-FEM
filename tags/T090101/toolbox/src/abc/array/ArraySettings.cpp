/* $Id: ArraySettings.cpp,v 1.1.1.1 2001-01-25 20:56:22 paklein Exp $ */
/* created: paklein (01/23/2001)                                          */
/* set the copy behavior for various templated array types                */

#include "ArrayT.h"

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
