/* $Id: ValueT.cpp,v 1.6 2003-04-22 18:32:16 paklein Exp $ */
#include "ValueT.h"
#include <stdlib.h>
#include <ctype.h>

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

/* enumeration */
ValueT::ValueT(const StringT& name, int value):
	fType(Enumeration),
	fInteger(value),
	fDouble(0.0)
{
	/* assign */
	operator=(name);
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
			out << fInteger;
			break;

		case Double:
			out << fDouble;
			break;

		case String:
			out << fString;
			break;

		case Enumeration:
		{
			/* write string if available */
			if (fString.StringLength() > 0)
				out << fString;
			else
				out << fInteger;
			break;
		}	
		default:
			ExceptionT::GeneralFail("ValueT::Write", "unsupported type %d", fType);
	}
}

namespace Tahoe {
ostream& operator<<(ostream& out, const ValueT& value)
{
	value.Write(out);
	return out;
}
}

ValueT& ValueT::operator=(int a)
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
	return *this;
}

ValueT& ValueT::operator=(double x)
{
	if (fType == Double)
		fDouble = x;
	else
		ExceptionT::GeneralFail("ValueT::operator=(double)", "type mismatch");	
	return *this;
}

ValueT& ValueT::operator=(const StringT& s)
{
	const char caller[] = "ValueT::operator=(StringT)";

	if (fType == String)
		fString = s;
	else if (fType == Enumeration)
	{
		fString = s;
		fString.Replace(' ', '_');

		/* checks */
		if (fString.StringLength() == 0)
			ExceptionT::GeneralFail(caller, "enumeration name cannot be empty");
		if (isdigit(fString[0]))
			ExceptionT::GeneralFail(caller, "enumeration name cannot start with [0-9]: \"%s\"",
				fString.Pointer());
	}
	else
		ExceptionT::GeneralFail(caller, "type mismatch");	
	return *this;
}

ValueT& ValueT::operator=(const ValueT& rhs)
{
	/* copy all values */
	fType = rhs.fType;
	fInteger = rhs.fInteger;
	fDouble = rhs.fDouble;
	fString = rhs.fString;
	
	/* dummy */
	return *this;
}

/* extract value from string, performing required type conversion */
void ValueT::FromString(const StringT& source)
{
	const char caller[] = "ValueT::FromString";

	/* cannot be empty */
	if (source.StringLength() == 0)
		ExceptionT::GeneralFail(caller, "source cannot be an empty string");

	switch (fType)
	{
		case Integer:
		{
			/* type conversion */
			int i = atoi(source.Pointer());
			
			/* assign */
			operator=(i);			
			break;
		}
		case Double:
		{
			/* type conversion */
			double d = atof(source.Pointer());
			
			/* assign */
			operator=(d);			
			break;
		}
		case String:
		{
			/* just copy */
			operator=(source);			
			break;
		}
		case Enumeration:
		{
			/* read number */
			if (isdigit(source[0])) 
			{	
				fType = Integer;
				FromString(source); /* treat as integer */
				fType = Enumeration;
			} 
			else 
			{
				fType = String;
				FromString(source); /* treat as string */
				fType = Enumeration;
			}
			break;
		}
		default:
			ExceptionT::GeneralFail(caller, "unsupported type %d", fType);
	}
}

/* type conversion operators not lvalues */
ValueT::operator const int&() const
{
	if (fType != Integer && fType != Enumeration)
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
	if (fType != String && fType != Enumeration)
		ExceptionT::GeneralFail("ValueT::operator const StringT&()", "type mismatch");	
	return fString;
}
