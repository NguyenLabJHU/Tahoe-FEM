/* $Id: LimitT.cpp,v 1.4 2003-04-22 18:32:16 paklein Exp $ */
#include "LimitT.h"

/* array behavior */
namespace Tahoe {
const bool ArrayT<LimitT>::fByteCopy = false;
}

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
