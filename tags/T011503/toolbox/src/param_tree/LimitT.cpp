/* $Id: LimitT.cpp,v 1.2 2002-09-03 07:54:08 paklein Exp $ */
#include "LimitT.h"

/* array behavior */
const bool ArrayT<LimitT>::fByteCopy = false;

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
