/* $Id: ValueT.cpp,v 1.5 2003-03-08 01:57:27 paklein Exp $ */
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

ValueT::ValueT(const ValueT& source):
	fType(source.fType),
	fInteger(source.fInteger),
	fDouble(source.fDouble),
	fString(source.fString)
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

namespace Tahoe {
ostream& operator<<(ostream& out, const ValueT& value)
{
	value.Write(out);
	return out;
}
}

int ValueT::operator=(int a)
{
	switch (fType)
	{
		case Integer:
		case Enumeration:
			fInteger = a;
			break;

		case Double:
			fDouble = double(a);
			break;

		default:
			ExceptionT::GeneralFail("ValueT::operator=(int)", "type mismatch");	
	}
	return a;
}

double ValueT::operator=(double x)
{
	if (fType == Double)
		fDouble = x;
	else
		ExceptionT::GeneralFail("ValueT::operator=(double)", "type mismatch");	
	return x;
}

const StringT& ValueT::operator=(const StringT& s)
{
	if (fType == String || fType == Enumeration)
		fString = s;
	else
		ExceptionT::GeneralFail("ValueT::operator=(StringT)", "type mismatch");	
	return s;
}

#if 0
const ValueT& ValueT::operator=(const ValueT& rhs)
{
	/* copy contents */
	fType = rhs.fType;
	fInteger = rhs.fInteger;
	fDouble = rhs.fDouble;
	fString = rhs.fString;
}
#endif

/* type conversion operators not lvalues */
ValueT::operator const int&() const
{
	if (fType != Integer)
		ExceptionT::GeneralFail("ValueT::operator const int&()", "type mismatch");	
	return fInteger;
}

ValueT::operator const double&() const
{
	if (fType != Double)
		ExceptionT::GeneralFail("ValueT::operator const double&()", "type mismatch");	
	return fDouble;
}

ValueT::operator const StringT&() const
{
	if (fType != String)
		ExceptionT::GeneralFail("ValueT::operator const StringT&()", "type mismatch");	
	return fString;
}
