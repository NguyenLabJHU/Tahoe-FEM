/* $Id: BandedLAdMatrixT.h,v 1.5 2003-11-21 22:41:36 paklein Exp $ */
/* created: MLK (05/21/1997)                                              */
/* square banded matrix operations                                        */
/* banded matrix elements stored in columns                               */

#ifndef _BANDED_LADMATRIX_T_H_
#define _BANDED_LADMATRIX_T_H_

/* Environmental */
#include "Environment.h"

/* base class */
#include "nArrayT.h"

namespace Tahoe {

/* forward declarations */
class dArrayT;
class dMatrixT;

class BandedLAdMatrixT: public nArrayT<double>
{
public:

	/** constructor */
	BandedLAdMatrixT(int squaredim, int leftbandsize, int rightbandsize);

	/** \name element and column accessor */
	/*@{*/
	double& operator()(int row, int col);
	const double& operator()(int row, int col) const;

	double* operator()(int col);
	const double* operator()(int col) const;

	/** returns 0.0 out of the band */
	double GetElement(int row, int col) const;
	/*@}*/

	/* assemble beginning with row and col in the upper left. */
	void AddBlock(int row, int col, const dMatrixT& block);
	
	/** \name dimensions */
	/*@{*/
	int Lband(void) const;
	int Rband(void) const;
	int Rows(void) const;
	int Cols(void) const;
	/*@}*/
	
	/* transpose copy */
	void Transpose(const BandedLAdMatrixT& matrix);

	/* I/O operator */
	friend ostream& operator<<(ostream& out, const BandedLAdMatrixT& matrix);

	/* Gaussian elimination with the given RHS vector or RHS matrix */
	void LinearSolve(dArrayT& RHS);
	void LinearSolve(dMatrixT& RHS);
	
	/*
	 * perform Gaussian elimination to compute the
	 * inverse of the banded matrix.  The result
	 * is stored in the RHS matrix.
	 */
	void BandedInverse(dMatrixT& RHS);
	
	/* bandedLAdMatrixT-vector multiplication - returns the result in b */
	void Multx(const dArrayT& x, dArrayT& b) const;
	void MultTx(const dArrayT& x, dArrayT& b) const;
	
	/* bandedLAdMatrixT-Matrix multiplication - returns the result in B */
	void MultM(const dMatrixT& M, dMatrixT& B) const;
	void MultTM(const dMatrixT& M, dMatrixT& B) const;
	
protected:

	int	fRows;
	int	fCols;
	int fLband;
	int fRband;
	int fColumnHeight;	
};

/* in-lines */

/* element accessor */
inline const double& BandedLAdMatrixT::operator()(int row, int col) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (row < 0 || row >= fRows || col < 0 || col >= fCols)
		throw ExceptionT::kOutOfRange;

	/* don't allow access to non banded elements */
	if(row-col > fLband || col-row > fRband) throw ExceptionT::kOutOfRange;
#endif
	
	return (fArray[col*fColumnHeight + fRband + row - col]);
};

inline double& BandedLAdMatrixT::operator()(int row, int col)
{
/* range checking */
#if __option (extended_errorcheck)
	if (row < 0 || row >= fRows || col < 0 || col >= fCols)
		throw ExceptionT::kOutOfRange;

	/* don't allow access to non banded elements */
	if(row-col > fLband || col-row > fRband) throw ExceptionT::kOutOfRange;
#endif
	
	return (fArray[col*fColumnHeight + fRband + row - col]);
};

/* returns a pointer to the top of the specified column */
inline double* BandedLAdMatrixT::operator()(int col)
{
/* range checking */
#if __option (extended_errorcheck)
	if (col < 0 || col >= fCols)
		throw ExceptionT::kOutOfRange;
#endif
	
	return (fArray + col*fColumnHeight);
};

const inline double* BandedLAdMatrixT::operator()(int col) const
{
/* range checking */
#if __option (extended_errorcheck)
	if (col < 0 || col >= fCols)
		throw ExceptionT::kOutOfRange;
#endif
	
	return (fArray + col*fColumnHeight);
};

/* dimensions */
inline int BandedLAdMatrixT::Lband(void) const { return (fLband); }
inline int BandedLAdMatrixT::Rband(void) const { return (fRband); }
inline int BandedLAdMatrixT::Rows(void) const { return (fRows); }
inline int BandedLAdMatrixT::Cols(void) const { return (fCols); }

} // namespace Tahoe 
#endif /* _BANDED_LADMATRIX_T_H_ */
