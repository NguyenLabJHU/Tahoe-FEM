/* $Id: LimitT.h,v 1.2 2002-09-03 07:54:08 paklein Exp $ */
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
		Lower,
		Upper,
		Only,
		Default
	};

	/** \name constructors */
	/*@{*/
	LimitT(int a, BoundT bound);
	LimitT(double x, BoundT bound);
	LimitT(const StringT& s, BoundT bound);
	/*@}*/
	
	/** return bound type */
	BoundT Bound(void) const { return fBound; };
	
	/** assess if the value satisfies the limit */
	bool InBound(const ValueT& value) const;

protected:

	/** bounde type */
	BoundT fBound;	
};

} // namespace Tahoe 
#endif /* _LIMIT_T_H_ */
