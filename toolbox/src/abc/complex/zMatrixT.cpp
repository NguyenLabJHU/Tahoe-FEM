/* $ Id $ */
/* created: paklein (05/19/1997)                                          */
/* 2 dimensional matrix mathematics object.                               */

#include "zMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "dMatrixT.h"

/*
* constructor
*/
zMatrixT::zMatrixT(void) { }
zMatrixT::zMatrixT(int numrows, int numcols): nMatrixT<ComplexT>(numrows,numcols) { }
zMatrixT::zMatrixT(int squaredim): nMatrixT<ComplexT>(squaredim) { }
zMatrixT::zMatrixT(int numrows, int numcols, ComplexT* p):
	nMatrixT<ComplexT>(numrows, numcols, p) { }
zMatrixT::zMatrixT(const dMatrixT& re, const dMatrixT& im)
{
	toZ(re,im);	
}
zMatrixT::zMatrixT(const zMatrixT& source):nMatrixT<ComplexT>(source) { }

/*
* I/O operator
*/
istream& operator>>(istream& in, zMatrixT& matrix)
{
	for (int j = 0; j < matrix.fRows; j++)
		for (int i = 0; i < matrix.fCols; i++)
				in >> matrix(j,i);

	return (in);
}

ostream& operator<<(ostream& out, const zMatrixT& matrix)
{
	for (int j = 0; j < matrix.fRows; j++)
	{
		for (int i = 0; i < matrix.fCols; i++)
				out << matrix(j,i);
		
		out << '\n';
	}
	
	return (out);
}

/*
* Returning the Real and Imaginary parts
*/
void zMatrixT::toRe(dMatrixT& re) const
{
	/* dimension check */
	if (fRows != re.Rows() || fCols != re.Cols()) throw(eGeneralFail);

	/* ComplexT function */
	z_to_Re(*this, re);
}

void zMatrixT::toIm(dMatrixT& im) const
{
	/* dimension check */
	if (fRows != im.Rows() || fCols != im.Cols()) throw(eGeneralFail);

	/* ComplexT function */
	z_to_Im(*this, im);
}

zMatrixT& zMatrixT::toZ(const dMatrixT& re, const dMatrixT& im)
{
	/* dimension checks */
	if (re.Rows() != re.Rows() ||
	    re.Cols() != im.Cols()) throw(eGeneralFail);
	
	/* dimension */
	Allocate(re.Rows(),im.Cols());
	
	/* ComplexT function */
	 ReIm_to_z(re,im,*this);

	return (*this);
}
