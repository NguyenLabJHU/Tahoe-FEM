/* $Id: LimitT.cpp,v 1.1 2002-09-03 07:04:33 paklein Exp $ */
#include "LimitT.h"

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

/* assess if the value satisfies the limit */
bool LimitT::InBound(const ValueT& value) const
{
#pragma message("LimitT::InBound: write me")
#pragma unused(value)
	return true;
}
