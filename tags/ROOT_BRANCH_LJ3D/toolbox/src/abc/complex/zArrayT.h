/* $Id: zArrayT.h,v 1.7 2003-11-21 22:41:33 paklein Exp $ */
/* created: PAK/AFLP (05/19/1997) */
#ifndef _ZARRAY_T_H_
#define _ZARRAY_T_H_

/* base class */
#include "nArrayT.h"

/* direct members */
#include "ComplexT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;

class zArrayT: public nArrayT<ComplexT>
{
public:

	/** constructors */
	/*@{*/
	zArrayT(void);
	explicit zArrayT(int length);
	zArrayT(int length, const ComplexT* p);
	zArrayT(const dArrayT& re, const dArrayT& im);
	zArrayT(const zArrayT& source);
	/*@}*/
	
	/*
	 * I/O operators
	 */
	friend istream& operator>>(istream& in, zArrayT& array);
	friend ostream& operator<<(ostream& out, const zArrayT& array);

	/*
	 * Assigment operators
	 */
	zArrayT& operator=(const zArrayT& RHS);
	zArrayT& operator=(const ComplexT& value);
	
	/*
	 * Returning the Real and Imaginary parts
	 */
	void toRe(dArrayT& re) const;
	void toIm(dArrayT& im) const;
	zArrayT& toZ(const dArrayT& re, const dArrayT& im);
	
	/* Conjugate every element in the array */
	zArrayT& Conjugate(const zArrayT& array);
	zArrayT& Conjugate(void);
};

/* Conjugate every element in this */
inline zArrayT& zArrayT::Conjugate(void) { return Conjugate(*this); }

/*
* Assigment operators
*/
inline zArrayT& zArrayT::operator=(const zArrayT& RHS)
{
	nArrayT<ComplexT>::operator=(RHS);
	return(*this);
}

inline zArrayT& zArrayT::operator=(const ComplexT& value)
{
	nArrayT<ComplexT>::operator=(value);
	return(*this);
}

} // namespace Tahoe 
#endif /* _ZARRAY_T_H_ */
