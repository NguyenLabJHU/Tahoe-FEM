/* $Id: zMatrixT.h,v 1.7 2003-11-21 22:41:33 paklein Exp $ */
/* created: paklein (05/19/1997) */
#ifndef _ZMATRIX_T_H_
#define _ZMATRIX_T_H_

/* base class */
#include "nMatrixT.h"

/* direct members */
#include "ComplexT.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;

/** two-dimensional matrix mathematics object */
class zMatrixT: public nMatrixT<ComplexT>
{
public:

	/** \name constructor */
	/*@{*/
	zMatrixT(void);
	zMatrixT(int numrows, int numcols);
	explicit zMatrixT(int squaredim);
	zMatrixT(int numrows, int numcols, const ComplexT* p);
	zMatrixT(const dMatrixT& re, const dMatrixT& im);
	zMatrixT(const zMatrixT& source);
	/*@}*/
	
	/*
	 * I/O operators
	 */
	friend istream& operator>>(istream& in, zMatrixT& matrix);
	friend ostream& operator<<(ostream& out, const zMatrixT& matrix);

	/*
	 * Assigment operators
	 */
	zMatrixT& operator=(const zMatrixT& RHS);
	zMatrixT& operator=(const ComplexT& value);

	/*
	 * Returning the Real and Imaginary parts
	 */
	void toRe(dMatrixT& re) const;
	void toIm(dMatrixT& im) const;
	zMatrixT& toZ(const dMatrixT& re, const dMatrixT& im);

  	/* matrix inverse */
  	zMatrixT& Inverse(const zMatrixT& matrix);
	zMatrixT& Inverse(void);

	/* conjugate every element in the matrix */
  	zMatrixT& Conjugate(const zMatrixT& matrix);  
	zMatrixT& Conjugate(void);  	
};

/* operate on self */
inline zMatrixT& zMatrixT::Inverse(void) { return Inverse(*this); }
inline zMatrixT& zMatrixT::Conjugate(void) { return Conjugate(*this); }

/*
* Assigment operators
*/
inline zMatrixT& zMatrixT::operator=(const zMatrixT& RHS)
{
	nMatrixT<ComplexT>::operator=(RHS);
	return(*this);
}

inline zMatrixT& zMatrixT::operator=(const ComplexT& value)
{
	nMatrixT<ComplexT>::operator=(value);
	return(*this);
}

} // namespace Tahoe 
#endif /* _ZMATRIX_T_H_ */
