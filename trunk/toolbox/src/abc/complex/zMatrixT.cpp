/* $Id: zMatrixT.cpp,v 1.9 2003-11-21 22:41:33 paklein Exp $ */
/* created: paklein (05/19/1997) */
#include "zMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "dMatrixT.h"

using namespace Tahoe;

zMatrixT::zMatrixT(void) { }
zMatrixT::zMatrixT(int numrows, int numcols): nMatrixT<ComplexT>(numrows,numcols) { }
zMatrixT::zMatrixT(int squaredim): nMatrixT<ComplexT>(squaredim) { }
zMatrixT::zMatrixT(int numrows, int numcols, const ComplexT* p):
	nMatrixT<ComplexT>(numrows, numcols, p) { }
zMatrixT::zMatrixT(const dMatrixT& re, const dMatrixT& im)
{
	toZ(re,im);	
}
zMatrixT::zMatrixT(const zMatrixT& source):nMatrixT<ComplexT>(source) { }

/*
* I/O operator
*/

namespace Tahoe {

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

} // namespace Tahoe

/*
* Returning the Real and Imaginary parts
*/
void zMatrixT::toRe(dMatrixT& re) const
{
	/* dimension check */
	if (fRows != re.Rows() || fCols != re.Cols()) throw ExceptionT::kOutOfRange;

	/* ComplexT function */
	ComplexT::z_to_Re(*this, re);
}

void zMatrixT::toIm(dMatrixT& im) const
{
	/* dimension check */
	if (fRows != im.Rows() || fCols != im.Cols()) throw ExceptionT::kOutOfRange;

	/* ComplexT function */
	ComplexT::z_to_Im(*this, im);
}

zMatrixT& zMatrixT::toZ(const dMatrixT& re, const dMatrixT& im)
{
	/* dimension checks */
	if (re.Rows() != re.Rows() ||
	    re.Cols() != im.Cols()) throw ExceptionT::kOutOfRange;
	
	/* dimension */
	Dimension(re.Rows(),im.Cols());
	
	/* ComplexT function */
	ComplexT::ReIm_to_z(re,im,*this);

	return (*this);
}


zMatrixT& zMatrixT::Inverse( const zMatrixT& matrix)
{
	/* dimension check */
	if (fRows != fCols || 
	   (fRows != 2 && fRows != 3)) throw ExceptionT::kSizeMismatch;
	
	/* (2 x 2) */
	if (fRows == 2)
	{
		/* temps - incase matrix is *this */
		ComplexT A0 = matrix.fArray[0];           
		ComplexT A1 = matrix.fArray[1];           
		ComplexT A2 = matrix.fArray[2];           
		ComplexT A3 = matrix.fArray[3];           

		ComplexT det = A0*A3 - 1*A1*A2;
		
				             
		fArray[0] =	A3/det;	             
		fArray[1] = -1*A1/det;
		fArray[2] = -1*A2/det;
		fArray[3] =	A0/det;		             
	}
	/* (3 x 3) */
	else
	{
		ComplexT z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13;
		ComplexT z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25;	

		z1 = matrix(0,0);
		z2 = matrix(0,1);
		z3 = matrix(0,2);
		z4 = matrix(1,0);
		z5 = matrix(1,1);
		z6 = matrix(1,2);
		z7 = matrix(2,0);
		z8 = matrix(2,1);
		z9 = matrix(2,2);
		z10 = -1.0*z2*z4;
		z11 = z3*z4;
		z12 = z1*z5;
		z13 = -1.0*z3*z5;
		z14 = -1.0*z1*z6;
		z15 = z2*z6;
		z16 = z13*z7;
		z17 = z15*z7;
		z18 = z2*z7;
		z19 = -1.0*z3*z7;
		z20 = -1.0*z5*z7;
		z7 = z6*z7;
		z21 = -1.0*z1*z8;
		z22 = z11*z8;
		z23 = z14*z8;
		z3 = z3*z8;
		z24 = z4*z8;
		z6 = -1.0*z6*z8;
		z1 = z1*z9;
		z8 = z10*z9;
		z25 = z12*z9;
		z2 = -1.0*z2*z9;
		z4 = -1.0*z4*z9;
		z5 = z5*z9;
		z9 = z10 + z12;
		z10 = z11 + z14;
		z11 = z13 + z15;
		z12 = z18 + z21;
		z13 = z20 + z24;
		z1 = z1 + z19;
		z8 = z16 + z17 + z22 + z23 + z25 + z8;
		z2 = z2 + z3;
		z3 = z4 + z7;
		z4 = z5 + z6;
		z5 = 1.0/z8;
		z6 = z5*z9;
		z7 = z10*z5;
		z8 = z11*z5;
		z9 = z12*z5;
		z10 = z13*z5;
		z1 = z1*z5;
		z2 = z2*z5;
		z3 = z3*z5;
		z4 = z4*z5;
		//{{z4, z2, z8},
		// {z3, z1, z7},
		// {z10, z9, z6}}

		ComplexT* pthis = Pointer();

		*pthis++ = z4;
		*pthis++ = z3;
		*pthis++ = z10;
		*pthis++ = z2;
		*pthis++ = z1;
		*pthis++ = z9;
		*pthis++ = z8;
		*pthis++ = z7;
		*pthis   = z6;
	}

	return(*this);
}

/* conjugate every element in the matrix */
zMatrixT& zMatrixT::Conjugate(const zMatrixT& matrix)
{
  /* must have same length */
  if (matrix.Length() != Length()) throw ExceptionT::kSizeMismatch;

  ComplexT* pLHS = Pointer();
  const ComplexT* pRHS = matrix.Pointer();
  int length = Length();
  for(int i = 0; i < length; i++)
	(*pLHS++).Conjugate(*pRHS++);
	
  return *this;
}
