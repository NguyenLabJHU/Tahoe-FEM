/* $Id: ValueT.cpp,v 1.3 2002-11-16 20:50:21 paklein Exp $ */
#include "ValueT.h"

/* array behavior */
namespace Tahoe {
const bool ArrayT<ValueT>::fByteCopy = false;
}

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

ValueT::ValueT(void):
	fType(None),
	fInteger(0),
	fDouble(0.0)
{

}

/* write the value to the output stream */
void ValueT::Write(ostream& out) const
{
	switch (fType)
	{
		case None:
			break;
			
		case Integer:
		case Enumeration:
			out << fInteger;
			break;

		case Double:
			out << fDouble;
			break;

		case String:
			out << fString;
			break;
	
		default:
			ExceptionT::GeneralFail();	
	}
}
