/* $Id: ParameterT.cpp,v 1.16 2005-04-05 15:51:45 paklein Exp $ */
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

ParameterT::ParameterT(const StringT& s, const char* name):
	ValueT(s.Pointer()),
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

	/* check type */
	if (fType != limit.Type()) {
		bool OK = false;
		switch (fType) {
			case Boolean:
				if (limit.Type() == Integer) OK = true;
				break;
			case Double:
				if (limit.Type() == Integer) OK = true;
				break;
			case String:
				if (limit.Type() == Word) OK = true;
				break;
			case Word:
				if (limit.Type() == String) OK = true;
				break;
		}
		if (!OK)
			ExceptionT::TypeMismatch(caller, "\"%s\" is of type \"%s\" not \"%s\"",
				Name().Pointer(), ValueT::TypeName(fType), ValueT::TypeName(limit.Type()));
	}

	bool added = fLimits.AppendUnique(limit);
	if (!added && fType == Enumeration) {
		StringT s = limit;
		if (s.StringLength() > 0)
			ExceptionT::GeneralFail(caller, "enumeration \"%s\" is not unique in \"%s\"", 
				s.Pointer(), Name().Pointer());
		else {
			int i = limit;
			ExceptionT::GeneralFail(caller, "enumeration %d is not unique in \"%s\"", 
				i, Name().Pointer());
		}
	}
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

	/* delete previous default */
	delete fDefault;

	switch (fType)
	{
		case Integer:
			fDefault = new ValueT(a);
			break;
		case Boolean:
			fDefault = new ValueT(bool(a));
			break;
		case Enumeration:
		{
			/* look for value in limits */
			for (int i = 0; i < fLimits.Length(); i++)
			{
				int i_limit = fLimits[i];
				if (i_limit == a) {
					fDefault = new ValueT(fLimits[i]);
					return;
				}
			}
			
			/* error on fall through */
			ExceptionT::GeneralFail(caller, "value %d does not appear in enumeration \"%s\"", a, fName.Pointer());
			break;
		}
		case Double:
			fDefault = new ValueT(double(a));
			break;
		default:
			ExceptionT::TypeMismatch(caller, "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
	}
}

void ParameterT::SetDefault(double x)
{
	/* delete previous default */
	delete fDefault;

	if (fType == Double)
		fDefault = new ValueT(x);
	else
		ExceptionT::TypeMismatch("ParameterT::SetDefault(double)", "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
}

void ParameterT::SetDefault(bool b)
{
	/* delete previous default */
	delete fDefault;

	if (fType == Boolean)
		fDefault = new ValueT(b);
	else
		ExceptionT::TypeMismatch("ParameterT::SetDefault(bool)", "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
}

void ParameterT::SetDefault(const char* s)
{
	const char caller[] = "ParameterT::SetDefault(StringT)";

	/* delete previous default */
	delete fDefault;

	/* check enumerations */
	if (fType == Enumeration)
	{
		/* look for value in limits */
		for (int i = 0; i < fLimits.Length(); i++)
		{
			const StringT& s_limit = fLimits[i];
			if (s_limit == s) 
			{
				fDefault = new ValueT(fLimits[i]);
				return;
			}
		}
			
		/* error if not found */
		ExceptionT::GeneralFail(caller, "value \"%s\" does not appear in enumeration \"%s\"", s, fName.Pointer());
	}
	else if (fType == String || fType == Word) {
		fDefault = new ValueT(s);

		/* check result of conversion to ValueT */
		if (fType != fDefault->Type())
			ExceptionT::GeneralFail(caller, "default \"%s\" is of type \"%s\" not \"%s\"",
				s, TypeName(fDefault->Type()), TypeName(fType));
	}
	else if (fType == Boolean)
	{	
		fDefault = new ValueT(Boolean);
		fDefault->FromString(s);
	}
	else
		ExceptionT::TypeMismatch(caller, "no conversion to %s for \"%s\"", TypeName(fType), fName.Pointer());
}

void ParameterT::SetDefault(const ValueT& v)
{
	switch (fType)
	{
		case Integer:
			SetDefault(int(v));
			break;
		case Double:
			SetDefault(double(v));
			break;
		case String:
		case Word:
		{
			const StringT& v_str = v;		
			SetDefault(v_str);
			break;
		}
		case Boolean:
			SetDefault(bool(v));
			break;
		case Enumeration:
		{
			/* try string, otherwise int */
			const StringT& v_str = v;
			if (v_str.StringLength() > 0)
				SetDefault(v_str);
			else
				SetDefault(int(v));
			break;
		}
		default:
			ExceptionT::GeneralFail("ParameterT::SetDefault");
	}
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
