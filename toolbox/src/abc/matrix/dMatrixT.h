/* $Id: dMatrixT.h,v 1.7 2002-07-02 19:56:47 cjkimme Exp $ */
/* created: paklein (05/24/1996) */

#ifndef _DMATRIX_T_H_
#define _DMATRIX_T_H_

/* base class */
#include "nMatrixT.h"


namespace Tahoe {

class dSymMatrixT;

/** double precision matrix. See nMatrixT for documentation on
 * inherited matrix operations and nArrayT for documentation on
 * additional mathematical operators */
class dMatrixT: public nMatrixT<double>
{
public:

	/** default constructor. Matrix has zero entries. */
	dMatrixT(void);

	/** constructor for generally dimensioned constructor. Matrix entries
	 * are not initialized. */
	dMatrixT(int numrows, int numcols);

	/** constructor for square matrix. Matrix entries are not initialized. 
	 * This constructor is declared explicit to avoid automatic int to
	 * dMatrixT type conversion */
	explicit dMatrixT(int squaredim);

	/** constructor for shallow matrix. Matrix will not free memory
	 * during destructor.
	 * \param numrows number of rows in the matrix
	 * \param numcols number of cols in the matrix
	 * \param p pointer to an array of length numrows*numcols */
	dMatrixT(int numrows, int numcols, double* p);

	/** copy constructor */
	dMatrixT(const dMatrixT& source);

	/* assignment operator. Operator will re-dimension matrix as needed.
	 * \param RHS source
	 * \return reference to *this */
	dMatrixT& operator=(const dMatrixT& RHS);

	/* assignment operator. Set all entries in the matrix to value
	 * \return reference to *this */
	dMatrixT& operator=(const double value);

	/** set this to the matrix inverse.
	 * \param matrix source matrix to inverse 
	 * \return reference to *this */
	dMatrixT& Inverse(const dMatrixT& matrix);

	/** set this its matrix inverse.
	 * \return reference to *this */
	dMatrixT& Inverse(void);

	/* matrix determinant.
	 * \note only implemented for (2 x 2) and (3 x 3) matrices */
	double Det(void) const;

	/* trace of the matrix.  
	 * \note The matrix must be square. */
	double Trace(void) const;

	/** returns the scalar inner product of a tensor with itself
	 * define as T:T = T_ij T_ij */
	double ScalarProduct(void) const { return Dot(*this, *this); };

	/* 2D/3D dimension transformations */
	void Rank2ExpandFrom2D(const dMatrixT& mat2D); /* fill with zeroes */
	void Rank2ReduceFrom3D(const dMatrixT& mat3D);

	/** set this to the symmetric part of matrix.
	 * \param matrix source matrix
	 * \return reference to *this */
	dMatrixT& Symmetrize(const dMatrixT& matrix);

	/** set this to its part.
	 * \return reference to *this */
	dMatrixT& Symmetrize(void);

	/** special matrix multiplication.
	 * \author thao */
	void MultSymAB(const dSymMatrixT& A, const dMatrixT& B);
	
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
}//namespace Tahoe
#endif /* _DMATRIX_T_H_ */
