/* $Id: ParameterT.h,v 1.7 2003-04-26 02:07:57 paklein Exp $ */
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
	ParameterT(int a, const StringT& name);
	ParameterT(double x, const StringT& name);
	ParameterT(const StringT& s, const StringT& name);

	/** set type without assigning value */
	ParameterT(TypeT t, const StringT& name);
	
	/** copy constructor */
	ParameterT(const ParameterT& source);
	
	/** default constructor. Needed to allow arrays of ParameterT's */
	ParameterT(void);
	/*@}*/
	
	/** destructor */
	~ParameterT(void);

	/** parameter name */
	const StringT& Name(void) const { return fName; };

	/** \name limits */
	/*@{*/
	/** add limit to parameter */
	void AddLimit(const LimitT& limit);

	void AddLimit(int a, LimitT::BoundT bound);
	void AddLimit(double x, LimitT::BoundT bound);
	void AddLimit(const StringT& s, LimitT::BoundT bound);

	/** return the list of limits */
	const ArrayT<LimitT>& Limits(void) const { return fLimits; };
	
	/** assess if the value satisties all limits */
	bool InBounds(const ValueT& value, bool verbose = false) const;
	/*@}*/

	/** \name set values with assignment operators 
	 * Only type conversion from int to double is allowed. All other
	 * type mismatched will through an exception. */
	/*@{*/
	ParameterT& operator=(int a);
	ParameterT& operator=(double x);
	ParameterT& operator=(const StringT& s);
	ParameterT& operator=(const ValueT& rhs);
	ParameterT& operator=(const ParameterT& rhs);
	/*@}*/

	/** \name description */
	/*@{*/
	void SetDescription(const StringT& description) { fDescription = description; };
	const StringT& Description(void) const { return fDescription; };
	/*@}*/

	/** \name default value */
	/*@{*/
	void SetDefault(int a);
	void SetDefault(double x);
	void SetDefault(const StringT& s);

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
inline void ParameterT::AddLimit(const StringT& s, LimitT::BoundT bound)
{
	LimitT limit(s, bound);
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
inline ParameterT& ParameterT::operator=(const StringT& s) { 
	ValueT::operator=(s); 
	return *this;
}
inline ParameterT& ParameterT::operator=(const ValueT& rhs) { 
	ValueT::operator=(rhs); 
	return *this;
}

} // namespace Tahoe 
#endif /* _PARAMETER_T_H_ */
