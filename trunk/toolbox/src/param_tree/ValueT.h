/* $Id: ValueT.h,v 1.4 2003-03-08 01:57:27 paklein Exp $ */
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
		Enumeration /**< stores values as integer or string */
	};

	/** \name constructors */
	/*@{*/
	ValueT(int a);
	ValueT(double x);
	ValueT(const StringT& s);

	/** set type without assigning value */
	ValueT(TypeT t);

	/** copy constructor */
	ValueT(const ValueT& source);
	
	/** default constructor */
	ValueT(void);
	/*@}*/

	/** write the value to the output stream */
	void Write(ostream& out) const;

	/** stream insertion operator */
	friend ostream& operator<<(ostream& out, const ValueT& value);

	/** value type */
	TypeT Type(void) const { return fType; };

	/** \name set values with assignment operators 
	 * Only type conversion from int to double is allowed. All other
	 * type mismatched will through an exception. */
	/*@{*/
	int operator=(int a);
	double operator=(double x);
	const StringT& operator=(const StringT& s);
	/*@}*/
	
	/** assignment operator */
//	const ValueT& operator=(const ValueT& rhs);
	
	/** \name type conversion operators not lvalues */
	/*@{*/
	operator const int&() const;
	operator const double&() const;
	operator const StringT&() const;
	/*@}*/

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
