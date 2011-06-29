/* $Id: ParameterT.h,v 1.12 2004-08-09 16:28:21 paklein Exp $ */
#ifndef _PARAMETER_T_H_
#define _PARAMETER_T_H_

/* base class */
#include "ValueT.h"

/* direct members */
#include "AutoArrayT.h"
#include "LimitT.h"

namespace Tahoe {

/** value with limit specifier */
class ParameterT: public ValueT
{
public:

	/** \name constructors */
	/*@{*/
	ParameterT(int a, const char* name);
	ParameterT(double x, const char* name);
	ParameterT(const char* s, const char* name);
	ParameterT(bool b, const char* name);
	ParameterT(const StringT& s, const char* name);

	/** set type without assigning value */
	ParameterT(TypeT t, const char* name);
	
	/** copy constructor */
	ParameterT(const ParameterT& source);
	
	/** default constructor. Needed to allow arrays of ParameterT's */
	ParameterT(void);
	/*@}*/
	
	/** destructor */
	~ParameterT(void);

	/** \name parameter name */
	/*@{*/
	const StringT& Name(void) const { return fName; };
	void SetName(const StringT& name) { fName = name; };
	/*@{*/

	/** \name limits */
	/*@{*/
	/** add limit to parameter */
	void AddLimit(const LimitT& limit);

	void AddLimit(int a, LimitT::BoundT bound);
	void AddLimit(double x, LimitT::BoundT bound);
	void AddLimit(const char* s, LimitT::BoundT bound);

	/** define a valid string-value pair */
	void AddEnumeration(const char* name, int value);

	/** add list of limits */
	void AddLimits(const ArrayT<LimitT>& limits);

	/** return the list of limits */
	const ArrayT<LimitT>& Limits(void) const { return fLimits; };
	
	/** assess if the value satisties all limits */
	bool InBounds(const ValueT& value, bool verbose = false) const;
	
	/** correct string-value pair. Fill in the missing string or value in
	 * the string-value pair in the given parameter using the enumerations 
	 * registered with ParameterT::AddEnumeration. */
	void FixEnumeration(ValueT& value) const;
	/*@}*/

	/** \name set values with assignment operators 
	 * Only type conversion from int to double is allowed. All other
	 * type mismatched will through an exception. */
	/*@{*/
	ParameterT& operator=(int a);
	ParameterT& operator=(double x);
	ParameterT& operator=(const char* s);
	ParameterT& operator=(const StringT& s);
	ParameterT& operator=(bool b);
	ParameterT& operator=(const ValueT& rhs);
	ParameterT& operator=(const ParameterT& rhs);
	/*@}*/

	/** \name description */
	/*@{*/
	void SetDescription(const char* description) { fDescription = description; };
	const StringT& Description(void) const { return fDescription; };
	/*@}*/

	/** \name default value */
	/*@{*/
	void SetDefault(int a);
	void SetDefault(double x);
	void SetDefault(const char* s);
	void SetDefault(const StringT& s);
	void SetDefault(bool b);
	void SetDefault(const ValueT& v);

	/** return a pointer to the default value or NULL if there isn't one */
	const ValueT* Default(void) const { return fDefault; };
	/*@}*/

protected:

	/** value name */
	StringT fName;
	
	/** description */
	StringT fDescription;

	/** default value */
	ValueT* fDefault;

	/** value limit specifications */
	AutoArrayT<LimitT> fLimits;
};

/* inlines */
inline void ParameterT::AddLimit(int a, LimitT::BoundT bound)
{
	LimitT limit(a, bound);
	AddLimit(limit);
}
inline void ParameterT::AddLimit(double x, LimitT::BoundT bound)
{
	LimitT limit(x, bound);
	AddLimit(limit);
}
inline void ParameterT::AddLimit(const char* s, LimitT::BoundT bound)
{
	LimitT limit(s, bound);
	AddLimit(limit);
}
inline void ParameterT::AddEnumeration(const char* name, int value)
{
	LimitT limit(name, value);
	AddLimit(limit);
}

inline ParameterT& ParameterT::operator=(int a) { 
	ValueT::operator=(a); 
	return *this;
}
inline ParameterT& ParameterT::operator=(double x) { 
	ValueT::operator=(x); 
	return *this;
}
inline ParameterT& ParameterT::operator=(const char* s) { 
	ValueT::operator=(s); 
	return *this;
}
inline ParameterT& ParameterT::operator=(bool b) { 
	ValueT::operator=(b); 
	return *this;
}
inline ParameterT& ParameterT::operator=(const ValueT& rhs) { 
	ValueT::operator=(rhs); 
	return *this;
}
inline ParameterT& ParameterT::operator=(const StringT& s) {
	return operator=(s.Pointer());
}

inline void ParameterT::SetDefault(const StringT& s) {
	SetDefault(s.Pointer());
}

} // namespace Tahoe 
#endif /* _PARAMETER_T_H_ */
