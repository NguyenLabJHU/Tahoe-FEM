/* $Id: ValueT.h,v 1.2 2002-11-16 20:50:21 paklein Exp $ */
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
		None,
		Integer,
		Double,
		String,
		Enumeration /**< stores values as integer */
	};

	/** \name constructors */
	/*@{*/
	ValueT(int a);
	ValueT(double x);
	ValueT(const StringT& s);

	/** set type without assigning value */
	ValueT(TypeT t);
	
	/** default constructor */
	ValueT(void);
	/*@}*/

	/** write the value to the output stream */
	void Write(ostream& out) const;

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
