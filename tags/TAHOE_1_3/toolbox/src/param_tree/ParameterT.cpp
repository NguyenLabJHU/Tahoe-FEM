/* $Id: ParameterT.cpp,v 1.11 2003-11-04 01:21:25 paklein Exp $ */
#include "ParameterT.h"

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ParameterT>::fByteCopy = false;
}

using namespace Tahoe;

/* constructors */
ParameterT::ParameterT(int a, const char* name):
	ValueT(a),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(double x, const char* name):
	ValueT(x),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(const char* s, const char* name):
	ValueT(s),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(bool b, const char* name):
	ValueT(b),
	fName(name),
	fDefault(NULL)
{

}

ParameterT::ParameterT(TypeT t, const char* name):
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
	const char caller[] = "ParameterT::AddLimit";

	/* check for enumerations */
	if (fType == Enumeration && limit.Bound() != LimitT::Only)
		ExceptionT::GeneralFail(caller, 
			"limits on enumerations for \"%s\" must be type \"only\"", fName.Pointer());

	/* no limits for booleans */
	if (fType == Boolean)
		ExceptionT::GeneralFail(caller, "no limits for boolean \"%s\"", fName.Pointer());	

	fLimits.Append(limit);
}

/* add list of limits */
void ParameterT::AddLimits(const ArrayT<LimitT>& limits)
{
	for (int i = 0; i < limits.Length(); i++)
		AddLimit(limits[i]);
}

/* correct string-value pair */
void ParameterT::FixEnumeration(ValueT& value) const
{
	const char caller[] = "ParameterT::FixEnumeration";
	if (fType != Enumeration || value.Type() != Enumeration)
		ExceptionT::TypeMismatch(caller);

	/* run through limits */
	for (int i = 0; i < fLimits.Length(); i++) {
		const LimitT& limit = fLimits[i];
		if (limit.InBound(value))
		{
			const StringT& value_s = value;

			/* fix integer value */ 
			if (value_s.StringLength() > 0) 
			{
				int limit_i = limit;
				value = limit_i;
			}
			/* fix string value */
			else 
			{
				const StringT& limit_s = limit;
				value = limit_s;
			}
			return;
		}
	}
			
	/* error on passing through */
	ExceptionT::GeneralFail(caller, "no matching value found");
}

/* assess if the value satisties all limits */
bool ParameterT::InBounds(const ValueT& value, bool verbose) const
{
	/* quick exit */
	if (fLimits.Length() == 0) return true;

	const char caller[] = "ParameterT::InBounds";
	
	/* flags for enumeration constraints */
	bool has_only = false;
	bool is_only = false;
	
	/* run through limits */
	for (int i = 0; i < fLimits.Length(); i++) {
	
		const LimitT& limit = fLimits[i];
		if (limit.Bound() == LimitT::Only) {
			has_only = true;
			if (!is_only) is_only = limit.InBound(value);
		}
		else if (!limit.InBound(value)) {

			if (verbose)
				cout << "\n " << caller << ": value " << value << " does not satisfy " 
					<< LimitT::ToString(limit.Bound()) << " bound " << limit << endl;

			return false;
		}
	}
	
	/* check enumeration limits */
	if (has_only) {
		
		/* message */
		if (!is_only && verbose)
			cout << "\n " << caller << ": value " << value << " does not satisfy " 
			     << LimitT::ToString(LimitT::Only) << " bounds" << endl;
	
		return is_only;	
	}
	else
		return true;
}

void ParameterT::SetDefault(int a)
{
	const char caller[] = "ParameterT::SetDefault(int)";

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
		case Boolean:
		{
			if (fDefault)
				*fDefault = bool(a);
			else
				fDefault = new ValueT(bool(a));
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
			ExceptionT::GeneralFail(caller, "value %d does not appear in enumeration \"%s\"", a, fName.Pointer());
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
			ExceptionT::TypeMismatch(caller, "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
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
		ExceptionT::TypeMismatch("ParameterT::SetDefault(double)", "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
}

void ParameterT::SetDefault(bool b)
{
	if (fType == Boolean)
	{
		if (fDefault)
			*fDefault = b;
		else
			fDefault = new ValueT(b);
	}
	else
		ExceptionT::TypeMismatch("ParameterT::SetDefault(bool)", "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
}

void ParameterT::SetDefault(const char* s)
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
			ExceptionT::GeneralFail(caller, "value \"%s\" does not appear in enumeration \"%s\"", s, fName.Pointer());
	}

	/* assign */
	if (fType == String || fType == Enumeration)
	{
		if (fDefault)
			*fDefault = s;
		else
			fDefault = new ValueT(s);
	}
	else if (fType == Boolean)
	{	
		if (!fDefault) fDefault = new ValueT(Boolean);
		fDefault->FromString(s);
	}
	else
		ExceptionT::TypeMismatch(caller, "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
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
