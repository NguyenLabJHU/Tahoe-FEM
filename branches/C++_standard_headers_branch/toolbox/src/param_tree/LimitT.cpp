/* $Id: LimitT.cpp,v 1.10 2004-09-14 18:13:31 paklein Exp $ */
#include "LimitT.h"

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<LimitT>::fByteCopy = false;

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

using namespace Tahoe;

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

LimitT::LimitT(const char* s, BoundT bound):
	ValueT(s),
	fBound(bound)
{

}

/* enumeration value */
LimitT::LimitT(const char* name, int value):
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
	switch (Type())
	{
		case Integer:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");
			
			int b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b > a;
			}
			else /* double */ {
				double a = value;
				return b > a;
			}
		}
		case Double:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");

			double b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b > a;
			}
			else /* double */ {
				double a = value;
				return b > a;
			}
		}
		case String:
		{
			/* type check */
			if (value.Type() != String)
				ExceptionT::GeneralFail(caller, "type mismatch");

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
	switch (Type())
	{
		case Integer:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");
			
			int b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b >= a;
			}
			else /* double */ {
				double a = value;
				return b >= a;
			}
		}
		case Double:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");

			double b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b >= a;
			}
			else /* double */ {
				double a = value;
				return b >= a;
			}
		}
		case String:
		{
			/* type check */
			if (value.Type() != String)
				ExceptionT::GeneralFail(caller, "type mismatch");

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
	switch (Type())
	{
		case Integer:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");
			
			int b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b < a;
			}
			else /* double */ {
				double a = value;
				return b < a;
			}
		}
		case Double:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");
			
			double b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b < a;
			}
			else /* double */ {
				double a = value;
				return b < a;
			}
		}
		case String:
		{
			/* type check */
			if (value.Type() != String)
				ExceptionT::GeneralFail(caller, "type mismatch");

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
	switch (Type())
	{
		case Integer:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");
			
			int b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b <= a;
			}
			else /* double */ {
				double a = value;
				return b <= a;
			}
		}
		case Double:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");

			double b = *this;
			if (value.Type() == Integer) {
				int a = value;
				return b <= a;
			}
			else /* double */ {
				double a = value;
				return b <= a;
			}
		}
		case String:
		{
			/* type check */
			if (value.Type() != String)
				ExceptionT::GeneralFail(caller, "type mismatch");

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
	switch (Type())
	{
		case Integer:
		{
			/* type check */
			if (value.Type() != Integer)
				ExceptionT::GeneralFail(caller, "type mismatch");

			int a = value;
			int b = *this;
			return a == b;
		}
		case Double:
		{
			/* type check */
			if (value.Type() != Integer && value.Type() != Double)
				ExceptionT::GeneralFail(caller, "type mismatch");

			double a = value;
			double b = *this;
			return a == b;
		}
		case String:
		{
			/* type check */
			if (value.Type() != String)
				ExceptionT::GeneralFail(caller, "type mismatch");

			const StringT& a = value;
			const StringT& b = *this;
			return a == b;
		}
		case Enumeration:
		{
			/* type check */
			if (value.Type() != Enumeration && value.Type() != String)
				ExceptionT::GeneralFail(caller, "type mismatch");

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
