/* $Id: dMatrixT.cpp,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (05/24/1996)                                          */

#include "dMatrixT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"

/* constructor */
dMatrixT::dMatrixT(void) { }
dMatrixT::dMatrixT(int numrows, int numcols): nMatrixT<double>(numrows,numcols) { }
dMatrixT::dMatrixT(int squaredim): nMatrixT<double>(squaredim) { }
dMatrixT::dMatrixT(int numrows, int numcols, double* p):
	nMatrixT<double>(numrows, numcols, p) { }
dMatrixT::dMatrixT(const dMatrixT& source): nMatrixT<double>(source) { }

/* matrix inverse functions - only implemented for (2 x 2)
* and (3 x 3) matrices */
dMatrixT& dMatrixT::Inverse(const dMatrixT& matrix)
{
	/* dimension check */
	if (fRows != fCols ||
	   (fRows != 2 && fRows != 3)) throw eSizeMismatch;
	
	/* (2 x 2) */
	if (fRows == 2)
	{
		/* temps - incase matrix is *this */
		double A0 = matrix.fArray[0];
		double A1 = matrix.fArray[1];
		double A2 = matrix.fArray[2];
		double A3 = matrix.fArray[3];

		double det = A0*A3 - A1*A2;
				
		fArray[0] =	A3/det;	
		fArray[1] =-A1/det;
		fArray[2] =-A2/det;
		fArray[3] =	A0/det;		
	}
	/* (3 x 3) */
	else
	{
		double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13;
		double z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24, z25;	

		z1 = matrix(0,0);
		z2 = matrix(0,1);
		z3 = matrix(0,2);
		z4 = matrix(1,0);
		z5 = matrix(1,1);
		z6 = matrix(1,2);
		z7 = matrix(2,0);
		z8 = matrix(2,1);
		z9 = matrix(2,2);
		z10 = -z2*z4;
		z11 = z3*z4;
		z12 = z1*z5;
		z13 = -z3*z5;
		z14 = -z1*z6;
		z15 = z2*z6;
		z16 = z13*z7;
		z17 = z15*z7;
		z18 = z2*z7;
		z19 = -z3*z7;
		z20 = -z5*z7;
		z7 = z6*z7;
		z21 = -z1*z8;
		z22 = z11*z8;
		z23 = z14*z8;
		z3 = z3*z8;
		z24 = z4*z8;
		z6 = -z6*z8;
		z1 = z1*z9;
		z8 = z10*z9;
		z25 = z12*z9;
		z2 = -z2*z9;
		z4 = -z4*z9;
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

		double* pthis = Pointer();

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

	return *this;
}

/* matrix determinants - only implemented for (2 x 2) and (3 x 3)
* matrices. */
double dMatrixT::Det(void) const
{
/* dimension check */
#if __option (extended_errorcheck)
	if (fRows != fCols ||
	   (fRows != 2 && fRows != 3)) throw eGeneralFail;
#endif
	
	if (fCols == 2) // (2 x 2)
		return fArray[0]*fArray[3] - fArray[1]*fArray[2];
	else            // (3 x 3)
		return fArray[0]*(fArray[4]*fArray[8] - fArray[5]*fArray[7])
			 - fArray[1]*(fArray[3]*fArray[8] - fArray[5]*fArray[6])
			 + fArray[2]*(fArray[3]*fArray[7] - fArray[4]*fArray[6]);
}

/* returns the Trace of the matrix.  Matrix must be square */
double dMatrixT::Trace(void) const
{
/* check is square */
#if __option (extended_errorcheck)
	if (fRows != fCols) throw eGeneralFail;
#endif

	double trace  = 0.0;

	if (fRows == 2)
	{
		trace += fArray[0];
		trace += fArray[3];
	}
	else if (fRows == 3)
	{
		trace += fArray[0];
		trace += fArray[4];
		trace += fArray[8];
	}
	else
	{
		double* pthis = fArray;
		int    offset = fRows + 1;
	
		for (int i = 0; i < fCols; i++)
		{
			trace += (*pthis);
			pthis += offset;
		}
	}
	
	return trace;
}

/* symmetrization */
dMatrixT& dMatrixT::Symmetrize(const dMatrixT& matrix)
{
#if __option (extended_errorcheck)	
	/* square matrices only */
	if (fRows != fCols ||
	    matrix.fRows != matrix.fCols ||
	    fRows != matrix.fRows) throw eSizeMismatch;
#endif

	if (fRows == 2)
		fArray[1] = fArray[2] = 0.5*(matrix.fArray[1] + matrix.fArray[2]);
	else if (fRows == 3)
	{
		fArray[1] = fArray[3] = 0.5*(matrix.fArray[1] + matrix.fArray[3]);
		fArray[2] = fArray[6] = 0.5*(matrix.fArray[2] + matrix.fArray[6]);
		fArray[5] = fArray[7] = 0.5*(matrix.fArray[5] + matrix.fArray[7]);	
	}
	else /* general routine */
	{
		for (int cols = 1; cols < fCols; cols++)
			for (int rows = 0; rows < cols; rows++)
				(*this)(rows,cols) = (*this)(cols,rows) =
					0.5*(matrix(rows,cols) + matrix(cols,rows));
	}
	
	return *this;
}

/* tranposition */
dMatrixT& dMatrixT::Transpose(const dMatrixT& matrix)
{
#if __option (extended_errorcheck)	
	if (fRows != matrix.fCols ||
	    fCols != matrix.fRows) throw eSizeMismatch;
#endif

	/* selve transposition */
	if (fArray == matrix.fArray) return Transpose();

	double *pthis = fArray;
	double *pm    = matrix.fArray;
	for (int i = 0; i < matrix.fRows; i++)
	{
		double* pmj = pm++;
	
		for (int j = 0; j < matrix.fCols; j++)
		{
			*pthis++ = *pmj;
			    pmj += matrix.fRows;
		}
	}
	
	return *this;
}

dMatrixT& dMatrixT::Transpose(void)
{
	register double temp;

	for (int i = 0; i < fRows - 1; i++)
	{
		double* prow = (*this)(i+1) + i;
		double* pcol = (*this)(i) + i + 1;
	
		for (int j = i + 1; j < fCols; j++)
		{
			temp  = *prow;
			*prow = *pcol;
			*pcol = temp;
			
			pcol++;
			prow += fRows;
		}
	}

	return *this;
}

/***********************************************
*
* Symmetric matrix specializations
*
**********************************************/

/* reduced index Rank 4 translations */
void dMatrixT::Rank4ReduceFrom3D(const dMatrixT& mat3D)
{
	/* dimension checks */
#if __option(extended_errorcheck)	
	if (fRows != fCols || fRows != 3) throw eGeneralFail;
	if (mat3D.fRows != mat3D.fCols || mat3D.fRows != 6) throw eSizeMismatch;
#endif

	double* pthis = fArray;
	
	*pthis++ = mat3D.fArray[0];	//1,1
	*pthis++ = mat3D.fArray[1]; //2,1
	*pthis++ = mat3D.fArray[5]; //3,1

	*pthis++ = mat3D.fArray[6]; //1,2
	*pthis++ = mat3D.fArray[7]; //2,2
	*pthis++ = mat3D.fArray[11]; //3,2

	*pthis++ = mat3D.fArray[30]; //3,1
	*pthis++ = mat3D.fArray[31]; //3,2
	*pthis   = mat3D.fArray[35]; //3,3
}

/* returns the Rank 4 devatoric operator in reduced index form and
* returns a reference to *this.
*
* Note: This operator cannot be used with a reduced index
*       vector to extract the deviatoric part by a simple Rank 2
*       matrix-vector operation because the terms in the vector
*       corresponding to the off-diagonal terms must be weighted
*       with a 2.  Use the dArrayT functions Deviatoric to do this */
dMatrixT& dMatrixT::ReducedIndexDeviatoric(void)
{
#if __option (extended_errorcheck)
	/* check */
	if (fRows != fCols || (fRows != 3 && fRows != 6)) throw eGeneralFail;
#endif

	*this = 0.0;
	
	double r23 = 2.0/3.0;
	double r13 =-1.0/3.0;

	if (fRows == 3)
	{
		(*this)(0,0) = (*this)(1,1) = r23;
		(*this)(0,1) = (*this)(1,0) = r13;
		(*this)(2,2) = 0.5;
	}
	else
	{
		(*this)(0,0) = (*this)(1,1) = (*this)(2,2) = r23;
		(*this)(0,1) = (*this)(1,0) = r13;
		(*this)(0,2) = (*this)(2,0) = r13;
		(*this)(2,1) = (*this)(1,2) = r13;
		(*this)(3,3) = (*this)(4,4) = (*this)(5,5) = 0.5;
	}

	return *this;
}

/* symmetric 4th rank tensor:
*
*	I_ijkl = 1/2 (d_ik d_jl + d_il d_jk)
*
* Returns a reference to this. */
dMatrixT& dMatrixT::ReducedIndexI(void)
{
#if __option (extended_errorcheck)
	/* check */
	if (fRows != fCols || (fRows != 3 && fRows != 6)) throw eGeneralFail;
#endif

	*this = 0.0;

	if (fRows == 3)
	{
		(*this)(0,0) = (*this)(1,1) = 1.0;
		(*this)(2,2) = 0.5;
	}
	else
	{
		(*this)(0,0) = (*this)(1,1) = (*this)(2,2) = 1.0;
		(*this)(3,3) = (*this)(4,4) = (*this)(5,5) = 0.5;
	}

	return *this;
}

/***********************************************
*
* Specializations added for element stiffness matrices - new class?
*
**********************************************/

/* expand into block diagonal submatrices if dimension factor */
void dMatrixT::Expand(const dMatrixT& B, int factor)
{
	/* dimension checks */
#if __option (extended_errorcheck)
	if (fRows != factor*B.fRows ||
	    fCols != factor*B.fCols) throw eSizeMismatch;
#endif

	double*	pCol  = Pointer();
	double* pBCol = B.Pointer();
	
	int coloffset = factor*fRows;
	
	for (int i = 0; i < B.fCols; i++)
	{
		double* pcol  = pCol;
		double* pBcol = pBCol;

		/* expand column of B */
		for (int k = 0; k < factor; k++)
		{
			double* psubBcol = pBcol;
			
			for (int j = 0; j < B.fRows; j++)
			{
				*pcol += *psubBcol++; /* accumulate */
				pcol  += factor;
			}
			
			pcol++;	/* shift down 1 in next column of this */
		}
			
		pCol  += coloffset;
		pBCol += B.fRows;
	}
}
