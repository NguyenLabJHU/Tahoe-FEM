/* $Id: ParameterT.cpp,v 1.1 2002-09-03 07:04:33 paklein Exp $ */
#include "ParameterT.h"

/* constructors */
ParameterT::ParameterT(int a, const StringT& name):
	ValueT(a),
	fName(name)
{

}

ParameterT::ParameterT(double x, const StringT& name):
	ValueT(x),
	fName(name)
{

}

ParameterT::ParameterT(const StringT& s, const StringT& name):
	ValueT(s),
	fName(name)
{

}

ParameterT::ParameterT(TypeT t, const StringT& name):
	ValueT(t),
	fName(name)
{

}
	
/* add limit to parameter */
void ParameterT::AddLimit(const LimitT& limit)
{
	fLimits.Append(limit);
}
