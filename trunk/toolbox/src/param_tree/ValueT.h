/* $Id: ValueT.h,v 1.6 2003-05-04 22:59:53 paklein Exp $ */
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
		None = 0,
		Integer,
		Double,
		String,
		Boolean,
		Enumeration /**< string-integer pair */
	};

	/** \name constructors */
	/*@{*/
	ValueT(int a);
	ValueT(double x);
	ValueT(const char* s);
	ValueT(bool b);

	/** enumeration. Enumerations are string-integer pairs. For all operators
	 * below, enumerations cast to both integers and strings. */
	ValueT(const char* s, int value);

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
	ValueT& operator=(int a);
	ValueT& operator=(double x);
	ValueT& operator=(bool b);
	ValueT& operator=(const char* s);
	ValueT& operator=(const StringT& s);
	ValueT& operator=(const ValueT& rhs);

	/** extract value from string, performing required type conversion */
	void FromString(const char* source);
	/*@}*/
	
	/** \name type conversion operators not lvalues */
	/*@{*/
	operator const int() const;
	operator const double() const;
	operator const bool() const;
	operator const StringT&() const;
	/*@}*/

	/** convert type name to string */
	static const char* TypeName(TypeT t);

protected:

	/** value type */
	TypeT fType;

	/** \name stored values */
	/*@{*/
	int fInteger;
	double fDouble;
	StringT fString;
	bool fBoolean;
	/*@}*/
};

/* inlines */
inline ValueT& ValueT::operator=(const StringT& s)
{
	return operator=(s.Pointer());
}

} /* namespace Tahoe */

#endif /* _VALUE_T_H_ */
