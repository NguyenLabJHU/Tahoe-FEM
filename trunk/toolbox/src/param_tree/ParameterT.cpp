/* $Id: ParameterT.cpp,v 1.4 2002-11-18 09:59:03 paklein Exp $ */
#include "ParameterT.h"

/* array behavior */
namespace Tahoe {
const bool ArrayT<ParameterT>::fByteCopy = false;
}

/* constructors */
ParameterT::ParameterT(int a, const StringT& name):
	ValueT(a),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(double x, const StringT& name):
	ValueT(x),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(const StringT& s, const StringT& name):
	ValueT(s),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(TypeT t, const StringT& name):
	ValueT(t),
	fName(name),
	fDefault(NULL)
{

}

/* copy constructor */
ParameterT::ParameterT(const ParameterT& source):
	ValueT(source),
	fDefault(NULL)
{
	/* duplicate default */
	if (source.fDefault) fDefault = new ValueT(*(source.fDefault));
}

/* default constructor */
ParameterT::ParameterT(void):fDefault(NULL) {}

/* destructor */
ParameterT::~ParameterT(void) { delete fDefault; }

/* add limit to parameter */
void ParameterT::AddLimit(const LimitT& limit)
{
	fLimits.Append(limit);
}

void ParameterT::SetDefault(int a)
{
	switch (fType)
	{
		case Integer:
		case Enumeration:
		{
			if (fDefault)
				*fDefault = a;
			else
				fDefault = new ValueT(a);
			break;
		}
		case Double:
		{
			if (fDefault)
				*fDefault = double(a);
			else
				fDefault = new ValueT(double(a));
			break;
		}
		default:
			ExceptionT::GeneralFail("ParameterT::SetDefault(int)", "type mismatch");	
	}
}

void ParameterT::SetDefault(double x)
{
	if (fType == Double)
	{
		if (fDefault)
			*fDefault = x;
		else
			fDefault = new ValueT(x);
	}
	else
		ExceptionT::GeneralFail("ParameterT::SetDefault(double)", "type mismatch");
}

void ParameterT::SetDefault(const StringT& s)
{
	if (fType == String || fType == Enumeration)
	{
		if (fDefault)
			*fDefault = s;
		else
			fDefault = new ValueT(s);
	}
	else
		ExceptionT::GeneralFail("ParameterT::SetDefault(StringT)", "type mismatch");
}

/* assignment operator */
const ParameterT& ParameterT::operator=(const ParameterT& rhs)
{
	/* inherited */
	ValueT::operator=(rhs);

	fName = rhs.fName;
	fDescription = rhs.fDescription;
	if (rhs.fDefault) {
		if (fDefault)
			*fDefault = *(rhs.fDefault);
		else
			fDefault = new ValueT(*(rhs.fDefault));
	}
	fLimits = rhs.fLimits;

	return *this;
}

