/* $Id: ArraySettings.cpp,v 1.2 2002-02-27 01:23:58 paklein Exp $ */
/* created: paklein (01/23/2001)                                          */
/* set the copy behavior for various templated array types                */

#include "ArrayT.h"

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
