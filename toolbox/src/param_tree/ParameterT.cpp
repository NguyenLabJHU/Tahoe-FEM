/* $Id: ParameterT.cpp,v 1.5 2003-04-22 18:32:16 paklein Exp $ */
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
	/* check for enumerations */
	if (fType == Enumeration && limit.Bound() != LimitT::Only)
		ExceptionT::GeneralFail("ParameterT::AddLimit", 
			"limits on enumerations must be type \"only\"");

	fLimits.Append(limit);
}

void ParameterT::SetDefault(int a)
{
	const char caller[] = "ParameterT::SetDefault";

	switch (fType)
	{
		case Integer:
		{
			if (fDefault)
				*fDefault = a;
			else
				fDefault = new ValueT(a);
			break;
		}
		case Enumeration:
		{
			/* look for value in limits */
			for (int i = 0; i < fLimits.Length(); i++)
			{
				int i_limit = fLimits[i];
				if (i_limit == a) {
					const StringT& s_limit = fLimits[i];
					SetDefault(s_limit);
					return;
				}
			}
			
			/* error on fall through */
			ExceptionT::GeneralFail(caller, "value %d does not appear in enumeration");
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
			ExceptionT::GeneralFail(caller, "type mismatch");	
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
	const char caller[] = "ParameterT::SetDefault(StringT)";

	/* check enumerations */
	if (fType == Enumeration)
	{
		/* look for value in limits */
		bool found = false;
		for (int i = 0; !found && i < fLimits.Length(); i++)
		{
			const StringT& s_limit = fLimits[i];
			found = (s_limit == s);
		}
			
		/* error */
		if (!found)
			ExceptionT::GeneralFail(caller, "value \"%s\" does not appear in enumeration",
				s.Pointer());
	}

	/* assign */
	if (fType == String || fType == Enumeration)
	{
		if (fDefault)
			*fDefault = s;
		else
			fDefault = new ValueT(s);
	}
	else
		ExceptionT::GeneralFail(caller, "type mismatch");
}

/* assignment operator */
ParameterT& ParameterT::operator=(const ParameterT& rhs)
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
