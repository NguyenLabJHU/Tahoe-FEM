/* $Id: ValueT.cpp,v 1.1 2002-09-03 07:04:33 paklein Exp $ */
#include "ValueT.h"

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
