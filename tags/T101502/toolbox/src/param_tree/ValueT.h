/* $Id: ValueT.h,v 1.1 2002-09-03 07:04:33 paklein Exp $ */
#ifndef _VALUE_T_H_
#define _VALUE_T_H_

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** basic parameter value */
class ValueT
{
public:

	/** enumerator for value type */
	enum TypeT {
		Integer,
		Double,
		String
	};

	/** \name constructors */
	/*@{*/
	ValueT(int a);
	ValueT(double x);
	ValueT(const StringT& s);

	/** set type without assigning value */
	ValueT(TypeT t);
	/*@}*/

	/** value type */
	TypeT Type(void) { return fType; };
	
protected:

	/** value type */
	TypeT fType;

	/** \name stored values */
	/*@{*/
	int fInteger;
	double fDouble;
	StringT fString;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _VALUE_T_H_ */
