/* $Id: ParameterT.h,v 1.3 2002-11-16 20:50:21 paklein Exp $ */
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
	
	/** default constructor. Should be used only to allow defining
	 * arrays of ParameterT's */
	ParameterT(void) {};
	/*@}*/

	/** parameter name */
	const StringT& Name(void) const { return fName; };

	/** \name limits */
	/*@{*/
	/** add limit to parameter */
	void AddLimit(const LimitT& limit);

	/** return the list of limits */
	const LinkedListT<LimitT>& Limits(void) const { return fLimits; } ;

	/** assess if the value satisties all limits */
	bool InBounds(const ValueT& value) const;
	/*@}*/

protected:

	/** value name */
	StringT fName;

	/** value limit specifications */
	LinkedListT<LimitT> fLimits;
};

} // namespace Tahoe 
#endif /* _PARAMETER_T_H_ */
