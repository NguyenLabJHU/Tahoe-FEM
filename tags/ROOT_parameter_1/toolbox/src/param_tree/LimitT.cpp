/* $Id: LimitT.cpp,v 1.5 2003-04-22 22:11:31 paklein Exp $ */
#include "LimitT.h"

/* array behavior */
namespace Tahoe {
const bool ArrayT<LimitT>::fByteCopy = false;

/* exceptions strings */
const char* LimitT::fBoundStrings[7] = 
{
/* 0 */ "no",
/* 1 */ "lower",
/* 2 */ "upper",
/* 3 */ "lower inclusive",
/* 4 */ "upper inclusive",
/* 5 */ "enumeration",
/* 6 */ "undefined"
};

/* return exception string */
const char* LimitT::ToString(BoundT bound)
{
	if (bound >= 0 && bound < 6)
		return fBoundStrings[bound];
	else
		return fBoundStrings[6];
}

} /* namespace Tahoe */

/* constructors */
LimitT::LimitT(int a, BoundT bound):
	ValueT(a),
	fBound(bound)
{

}

LimitT::LimitT(double x, BoundT bound):
	ValueT(x),
	fBound(bound)
{

}

LimitT::LimitT(const StringT& s, BoundT bound):
	ValueT(s),
	fBound(bound)
{

}

/* enumeration value */
LimitT::LimitT(const StringT& name, int value):
	ValueT(name, value),
	fBound(Only)
{

}

/* assess if the value satisfies the limit */
bool LimitT::InBound(const ValueT& value) const
{
	const char caller[] = "LimitT::InBound";

	/* bound type */
	switch (fBound)
	{
		case None:
			return true;

		case Lower:
			return CheckLower(value);

		case Upper:
			return CheckUpper(value);

		case LowerInclusive:
			return CheckLowerInclusive(value);

		case UpperInclusive:
			return CheckUpperInclusive(value);

		case Only:
			return CheckOnly(value);
		
		default:
			ExceptionT::GeneralFail(caller, "unrecognized bound");
	}

	/* catch all */
	return false;
}

/**********************************************************************
 * Private
 **********************************************************************/

bool LimitT::CheckUpper(const ValueT& value) const
{
	const char caller[] = "LimitT::CheckUpper";
	if (value.Type() != Type()) ExceptionT::GeneralFail(caller, "type mismatch");

	switch (Type())
	{
		case Integer:
		{
			int a = value;
			int b = *this;
			return b > a;
		}
		case Double:
		{
			double a = value;
			double b = *this;
			return b > a;
		}
		case String:
		{
			const StringT& a = value;
			const StringT& b = *this;
			return b > a;
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized type");
	}
	
	return false;
}

bool LimitT::CheckUpperInclusive(const ValueT& value) const
{
	const char caller[] = "LimitT::CheckUpperInclusive";
	if (value.Type() != Type()) ExceptionT::GeneralFail(caller, "type mismatch");

	switch (Type())
	{
		case Integer:
		{
			int a = value;
			int b = *this;
			return b >= a;
		}
		case Double:
		{
			double a = value;
			double b = *this;
			return b >= a;
		}
		case String:
		{
			const StringT& a = value;
			const StringT& b = *this;
			return (a == b || b > a);
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized type");
	}
	
	return false;
}

bool LimitT::CheckLower(const ValueT& value) const
{
	const char caller[] = "LimitT::CheckLower";
	if (value.Type() != Type()) ExceptionT::GeneralFail(caller, "type mismatch");

	switch (Type())
	{
		case Integer:
		{
			int a = value;
			int b = *this;
			return b < a;
		}
		case Double:
		{
			double a = value;
			double b = *this;
			return b < a;
		}
		case String:
		{
			const StringT& a = value;
			const StringT& b = *this;
			return b < a;
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized type");
	}
	
	return false;
}

bool LimitT::CheckLowerInclusive(const ValueT& value) const
{
	const char caller[] = "LimitT::CheckLowerInclusive";
	if (value.Type() != Type()) ExceptionT::GeneralFail(caller, "type mismatch");

	switch (Type())
	{
		case Integer:
		{
			int a = value;
			int b = *this;
			return b <= a;
		}
		case Double:
		{
			double a = value;
			double b = *this;
			return b <= a;
		}
		case String:
		{
			const StringT& a = value;
			const StringT& b = *this;
			return (a == b || b < a);
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized type");
	}
	
	return false;
}

bool LimitT::CheckOnly(const ValueT& value) const
{
	const char caller[] = "LimitT::CheckOnly";
	if (value.Type() != Type()) ExceptionT::GeneralFail(caller, "type mismatch");

	switch (Type())
	{
		case Integer:
		{
			int a = value;
			int b = *this;
			return a == b;
		}
		case Double:
		{
			double a = value;
			double b = *this;
			return a == b;
		}
		case String:
		{
			const StringT& a = value;
			const StringT& b = *this;
			return a == b;
		}
		case Enumeration:
		{
			const StringT& sa = value;
			
			/* try to compare string first */
			if (sa.StringLength() > 0)
			{
				const StringT& sb = *this;
				return sa == sb;
			}
			else /* compare integers */
			{
				int ia = value;
				int ib = *this;
				return ia == ib;
			}
		}
		default:
			ExceptionT::GeneralFail(caller, "unrecognized type");
	}
	
	return false;
}
