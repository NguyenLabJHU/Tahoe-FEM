/* $Id: LimitT.cpp,v 1.3 2003-03-14 23:45:09 paklein Exp $ */
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

/* assess if the value satisfies the limit */
bool LimitT::InBound(const ValueT& value) const
{
#pragma message("LimitT::InBound: write me")
#pragma unused(value)
	return true;
}
