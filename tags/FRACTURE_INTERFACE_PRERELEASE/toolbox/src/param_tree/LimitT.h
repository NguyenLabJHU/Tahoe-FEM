/* $Id: LimitT.h,v 1.6 2003-05-04 22:59:53 paklein Exp $ */
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
		None = 0,   /**< non-limit, needed for default constructor */
		Lower = 1,
		Upper = 2,
		LowerInclusive = 3,
		UpperInclusive = 4,
		Only = 5 /**< for fixed number of allowed values or enumerations */
	};

	/** \name constructors */
	/*@{*/
	LimitT(int a, BoundT bound);
	LimitT(double x, BoundT bound);
	LimitT(const char* s, BoundT bound);
	
	/** enumeration value. String-integer pair that bounds as LimitT::Only */
	LimitT(const char* name, int value);
	
	/** default constructor */
	LimitT(void): fBound(None) {};
	/*@}*/
	
	/** return bound type */
	BoundT Bound(void) const { return fBound; };
	
	/** assess if the value satisfies the limit */
	bool InBound(const ValueT& value) const;

	/** string associated with each bound type */
	static const char* ToString(BoundT bound);

private:

	/** \name bounds tests
	 * Return true if value satisfies bound */
	/*@{*/
	bool CheckLower(const ValueT& value) const;
	bool CheckLowerInclusive(const ValueT& value) const;
	bool CheckUpper(const ValueT& value) const;
	bool CheckUpperInclusive(const ValueT& value) const;
	bool CheckOnly(const ValueT& value) const;
	/*@}*/

protected:

	/** bound type */
	BoundT fBound;

  	/** bound types */
  	static const char* fBoundStrings[7];
};

} // namespace Tahoe 
#endif /* _LIMIT_T_H_ */
