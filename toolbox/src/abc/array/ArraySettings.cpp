/* $Id: ArraySettings.cpp,v 1.4 2002-07-02 19:56:39 cjkimme Exp $ */
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
