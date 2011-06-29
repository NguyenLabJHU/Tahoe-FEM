/* $Id: ValueT.cpp,v 1.14 2004-07-20 06:12:38 paklein Exp $ */
#include "ValueT.h"
#include <stdlib.h>
#include <ctype.h>

/* exceptions strings */
static const char* type_names[7] = {
/* 0 */ "none",
/* 1 */ "integer",
/* 2 */ "double",
/* 3 */ "string",
/* 4 */ "boolean",
/* 5 */ "enumeration",
/* 6 */ "word"};

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ValueT>::fByteCopy = false;
}

using namespace Tahoe;

/* constructors */
ValueT::ValueT(int a): 
	fType(Integer),
	fInteger(a),
	fDouble(0.0),
	fBoolean(false)
{

}

ValueT::ValueT(double x):
	fType(Double),
	fInteger(0),
	fDouble(x),
	fBoolean(false)
{

}

ValueT::ValueT(const char* s):
	fInteger(0),
	fDouble(0.0),
	fString(s),
	fBoolean(false)
{
	/* look for white space */
	while (*s != '\0') {
		if (isspace(*s)) {
			fType = String;
			return;
		}
		s++;
	}

	/* no space */
	fType = Word;
}

ValueT::ValueT(bool b):
	fType(Boolean),
	fInteger(0),
	fDouble(0.0),
	fBoolean(b)
{

}

/* enumeration */
ValueT::ValueT(const char* name, int value):
	fType(Enumeration),
	fInteger(value),
	fDouble(0.0),
	fBoolean(false)
{
	/* assign */
	operator=(name);
}

ValueT::ValueT(TypeT t):
	fType(t),
	fInteger(0),
	fDouble(0.0),
	fBoolean(false)
{

}

ValueT::ValueT(const ValueT& source):
	fType(source.fType),
	fInteger(source.fInteger),
	fDouble(source.fDouble),
	fString(source.fString),
	fBoolean(source.fBoolean)
{

}

ValueT::ValueT(void):
	fType(None),
	fInteger(0),
	fDouble(0.0),
	fBoolean(false)
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
		case Word:
			out << fString;
			break;

		case Boolean:
			out << ((fBoolean) ? "true" : "false");
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

		case Boolean:
			fBoolean = bool(a);
			break;

		case Double:
			fDouble = double(a);
			break;

		default:
			ExceptionT::TypeMismatch("ValueT::operator=(int)");	
	}
	return *this;
}

ValueT& ValueT::operator=(double x)
{
	if (fType == Double)
		fDouble = x;
	else
		ExceptionT::TypeMismatch("ValueT::operator=(double)");	
	return *this;
}

ValueT& ValueT::operator=(bool b)
{
	if (fType == Boolean)
		fBoolean = b;
	else
		ExceptionT::TypeMismatch("ValueT::operator=(bool)");	
	return *this;
}

ValueT& ValueT::operator=(const char* s)
{
	const char caller[] = "ValueT::operator=(const char*)";

	if (fType == String)
		fString = s;
	else if (fType == Word) {
		int count;
		fString.FirstWord(s, count, false);
	}
	else if (fType == Boolean)
		FromString(s);
	else if (fType == Enumeration)
	{
		fString = s;
		fString.Replace(' ', '_');

		/* checks */
		if (fString.StringLength() == 0)
			ExceptionT::GeneralFail(caller, "enumeration name cannot be empty");
#if 0
//TEMP - why is this name invalid?
		if (isdigit(fString[0]))
			ExceptionT::GeneralFail(caller, "enumeration name cannot start with [0-9]: \"%s\"",
				fString.Pointer());
#endif
	}
	else
		ExceptionT::TypeMismatch(caller);	
	return *this;
}

ValueT& ValueT::operator=(const ValueT& rhs)
{
	/* copy all values */
	fType = rhs.fType;
	fInteger = rhs.fInteger;
	fDouble = rhs.fDouble;
	fString = rhs.fString;
	fBoolean = rhs.fBoolean;
	
	/* dummy */
	return *this;
}

/* extract value from string, performing required type conversion */
void ValueT::FromString(const char* source)
{
	const char caller[] = "ValueT::FromString";

	/* cannot be empty */
	if (strlen(source) == 0)
		ExceptionT::GeneralFail(caller, "source cannot be an empty string");

	switch (fType)
	{
		case Integer:
		{
			/* type conversion */
			int i = atoi(source);
			
			/* assign */
			operator=(i);			
			break;
		}
		case Double:
		{
			/* type conversion */
			double d = atof(source);
			
			/* assign */
			operator=(d);			
			break;
		}
		case String:
		case Word:
		{
			/* just copy */
			operator=(source);			
			break;
		}
		case Boolean:
		{
			if (source[0] == 't' || source[0] == 'T' || source[0] == '1')
				fBoolean = true;
			else if (source[0] == 'f' || source[0] == 'F' || source[0] == '0')
				fBoolean = false;
			else
				ExceptionT::GeneralFail(caller, "could not extract bool from \"%s\"", source);
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
				fType = Word;
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
ValueT::operator const int() const
{
	if (fType == Integer || fType == Enumeration)
		return fInteger;
	else if (fType == Double)
		return int(fDouble);
	else if (fType == Boolean)
		return int(fBoolean);
	else
		ExceptionT::TypeMismatch("ValueT::int()", "cannot convert from %s", TypeName(fType));
		
	return 0;
}

ValueT::operator const bool() const
{
	if (fType != Boolean)
		ExceptionT::TypeMismatch("ValueT::bool()", "cannot convert from %s", TypeName(fType));
	return fBoolean;
}

ValueT::operator const double() const
{
	if (fType == Double)
		return fDouble;
	else if (fType == Integer)
		return double(fInteger);
	else
		ExceptionT::TypeMismatch("ValueT::double()", "cannot convert from %s", TypeName(fType));

	return 0.0;
}

ValueT::operator const StringT&() const
{
	if (fType != String && fType != Word && fType != Enumeration)
		ExceptionT::TypeMismatch("ValueT::StringT&()", "cannot convert from %s", TypeName(fType));
	return fString;
}

/* convert type name to string */
const char* ValueT::TypeName(TypeT t)
{
	if (t >= None && t <= Word)
		return type_names[t];
	else
		return type_names[0];
}

/* comparison */
bool ValueT::operator==(const ValueT& rhs) const
{
	switch (fType)
	{
		case Integer:
		{
			int i = rhs;
			return fInteger == i;
		}
		case Double:
		{
			double a = rhs;
			return fDouble == a;
		}
		case String:
		case Word:
		{
			const StringT& s = rhs;
			return fString == s;
		}
		case Boolean:
		{
			bool t = rhs;
			return fBoolean == t;
		}
		case Enumeration:
		{
			const StringT& s = rhs;
			if (s.StringLength() > 0)
				return fString == s;
			else {
				int i = rhs;
				return fInteger == i;
			}		
		}	
		default:
			return false;
	}
}
