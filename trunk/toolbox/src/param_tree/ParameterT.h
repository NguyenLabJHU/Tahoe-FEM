/* $Id: ParameterT.h,v 1.4 2002-11-18 09:59:03 paklein Exp $ */
#ifndef _PARAMETER_T_H_
#define _PARAMETER_T_H_

/* base class */
#include "ValueT.h"

/* direct members */
#include "LinkedListT.h"
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

	/** return the list of limits */
	const LinkedListT<LimitT>& Limits(void) const { return fLimits; };
	
	/** return non-const version of the list. Needed to traverse the list */
	LinkedListT<LimitT>& Limits(void) { return fLimits; };

	/** assess if the value satisties all limits */
	bool InBounds(const ValueT& value) const;
	/*@}*/

	/** \name set values with assignment operators 
	 * Only type conversion from int to double is allowed. All other
	 * type mismatched will through an exception. */
	/*@{*/
	int operator=(int a);
	double operator=(double x);
	const StringT& operator=(const StringT& s);
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

	/** assignment operator */
	const ParameterT& operator=(const ParameterT& rhs);

protected:

	/** value name */
	StringT fName;
	
	/** description */
	StringT fDescription;

	/** default value */
	ValueT* fDefault;

	/** value limit specifications */
	LinkedListT<LimitT> fLimits;
};

/* inlines */
inline int ParameterT::operator=(int a) { return ValueT::operator=(a); }
inline double ParameterT::operator=(double x) { return ValueT::operator=(x); }
inline const StringT& ParameterT::operator=(const StringT& s) { return ValueT::operator=(s); }

} // namespace Tahoe 
#endif /* _PARAMETER_T_H_ */
