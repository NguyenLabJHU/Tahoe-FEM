/* $Id: dMatrixT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (05/24/1996)                                          */

#ifndef _DMATRIX_T_H_
#define _DMATRIX_T_H_

/* base class */
#include "nMatrixT.h"

class dMatrixT: public nMatrixT<double>
{
public:

	/* constructor */
	dMatrixT(void);
	dMatrixT(int numrows, int numcols);
	dMatrixT(int squaredim);
	dMatrixT(int numrows, int numcols, double* p);
	dMatrixT(const dMatrixT& source);

	/* assignment operators */
	dMatrixT& operator=(const dMatrixT& RHS);
	dMatrixT& operator=(const double value);

	/* matrix inverse functions - only implemented for (2 x 2)
	 * and (3 x 3) matrices */
	dMatrixT& Inverse(const dMatrixT& matrix);
	dMatrixT& Inverse(void);

	/* matrix determinants - only implemented for (2 x 2) and (3 x 3)
	 * matrices */
	double Det(void) const;

	/* returns the Trace of the matrix.  Matrix must be square */
	double Trace(void) const;

	/* 2D/3D dimension transformations */
	void Rank2ExpandFrom2D(const dMatrixT& mat2D); /* fill with zeroes */
	void Rank2ReduceFrom3D(const dMatrixT& mat3D);

	/* symmetrization */
	dMatrixT& Symmetrize(const dMatrixT& matrix);
	dMatrixT& Symmetrize(void);

	/* transposition */
	dMatrixT& Transpose(const dMatrixT& matrix);
	dMatrixT& Transpose(void);
	
/***********************************************
* Symmetric matrix specializations
**********************************************/
	
	/* reduced index Rank 4 translations */
	void Rank4ReduceFrom3D(const dMatrixT& mat3D);

	/* returns the Rank 4 devatoric operator in reduced index form and
	 * returns a reference to *this.
	 *
	 * Note: This operator cannot be used with a reduced index
	 *       vector to extract the deviatoric part by a simple Rank 2
	 *       matrix-vector operation because the terms in the vector
	 *       corresponding to the off-diagonal terms must be weighted
	 *       with a 2.  Use the dArrayT functions Deviatoric to do this */
	dMatrixT& ReducedIndexDeviatoric(void);

	/*
	 * Symmetric 4th rank tensor:
	 *
	 *	I_ijkl = 1/2 (d_ik d_jl + d_il d_jk)
	 *
	 * Returns a reference to this.
	 *
	 */
	dMatrixT& ReducedIndexI(void);

/***********************************************
* Specializations added for element stiffness matrices - new class?
**********************************************/
	
	/* expand into block by block diagonal submatrices - factor made
	 * argument to save division to calculate it */
	void Expand(const dMatrixT& B, int factor); 	
};

/* inlines */

/* assigment operators */
inline dMatrixT& dMatrixT::operator=(const dMatrixT& RHS)
{
	nMatrixT<double>::operator=(RHS);
	return *this;
}

inline dMatrixT& dMatrixT::operator=(const double value)
{
	nMatrixT<double>::operator=(value);
	return *this;
}

inline dMatrixT& dMatrixT::Inverse(void) { return Inverse(*this); }
inline dMatrixT& dMatrixT::Symmetrize(void) { return Symmetrize(*this); }

/* 2D/3D dimension transformations */
inline void dMatrixT::Rank2ExpandFrom2D(const dMatrixT& mat2D)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != fCols || fRows != 3) throw eGeneralFail;
	if (mat2D.fRows != mat2D.fCols || mat2D.fRows != 2) throw eSizeMismatch;
#endif

	double* p = fArray;
	*p++ = mat2D.fArray[0];				
	*p++ = mat2D.fArray[1];
	*p++ = 0.0;

	*p++ = mat2D.fArray[2];
	*p++ = mat2D.fArray[3];
	*p++ = 0.0;

	*p++ = 0.0;
	*p++ = 0.0;
	*p   = 0.0;
}

inline void dMatrixT::Rank2ReduceFrom3D(const dMatrixT& mat3D)
{
	/* dimension checks */
#if __option (extended_errorcheck)	
	if (fRows != fCols || fRows != 2) throw eGeneralFail;
	if (mat3D.fRows != mat3D.fCols || mat3D.fRows != 3) throw eSizeMismatch;
#endif

	double* p = fArray;
	*p++ = mat3D.fArray[0];				
	*p++ = mat3D.fArray[1];

	*p++ = mat3D.fArray[3];
	*p   = mat3D.fArray[4];
}

#endif /* _DMATRIX_T_H_ */
