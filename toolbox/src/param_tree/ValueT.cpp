/* $Id: ValueT.cpp,v 1.2 2002-09-03 07:54:08 paklein Exp $ */
#include "ValueT.h"

/* array behavior */
const bool ArrayT<ValueT>::fByteCopy = false;

/* constructors */
ValueT::ValueT(int a): 
	fType(Integer),
	fInteger(a),
	fDouble(0.0)
{

}

ValueT::ValueT(double x):
	fType(Double),
	fInteger(0),
	fDouble(x)
{

}

ValueT::ValueT(const StringT& s):
	fType(String),
	fInteger(0),
	fDouble(0.0),
	fString(s)
{

}

ValueT::ValueT(TypeT t):
	fType(t),
	fInteger(0),
	fDouble(0.0)
{

}
