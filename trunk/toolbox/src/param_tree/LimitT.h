/* $Id: LimitT.h,v 1.4 2003-04-22 18:32:16 paklein Exp $ */
#ifndef _LIMIT_T_H_
#define _LIMIT_T_H_

/* base class */
#include "ValueT.h"

namespace Tahoe {

/** basic parameter value */
class LimitT: public ValueT
{
public:

	/** enumerator for limit type */
	enum BoundT {
	    None,   /**< non-limit, needed for default constructor */
		Lower,
		Upper,
		Only,   /**< for fixed number of allowed values or enumerations */
		Default /**< needed? ParameterT has separate field for default */
	};

	/** \name constructors */
	/*@{*/
	LimitT(int a, BoundT bound);
	LimitT(double x, BoundT bound);
	LimitT(const StringT& s, BoundT bound);
	
	/** enumeration value. String-integer pair that bounds as LimitT::Only */
	LimitT(const StringT& name, int value);
	
	/** default constructor */
	LimitT(void): fBound(None) {};
	/*@}*/
	
	/** return bound type */
	BoundT Bound(void) const { return fBound; };
	
	/** assess if the value satisfies the limit */
	bool InBound(const ValueT& value) const;

private:

	/** \name bounds tests
	 * Return true if value satisfies bound */
	/*@{*/
	bool CheckLower(const ValueT& value) const;
	bool CheckUpper(const ValueT& value) const;
	bool CheckOnly(const ValueT& value) const;
	/*@}*/

protected:

	/** bounde type */
	BoundT fBound;	
};

} // namespace Tahoe 
#endif /* _LIMIT_T_H_ */
