/* $Id: LAdMatrixT.h,v 1.1.1.1 2001-01-25 20:56:23 paklein Exp $ */
/* created: paklein (12/05/1996)                                          */
/* dMatrixT with some linear algebra functions                            */

#ifndef _LA_DMATRIX_T_H_
#define _LA_DMATRIX_T_H_

/* base class */
#include "dMatrixT.h"

/* forward declarations */
class dArrayT;

class LAdMatrixT: public dMatrixT
{
public:

	/* constructor */
	LAdMatrixT(void);
	LAdMatrixT(int squaredim);
	LAdMatrixT(const LAdMatrixT& source);

	/* pivoting functions */
	void RowPivot(int row1, int row2);
	void ColumnPivot(int col1, int col2);
	void SymmetricPivot(int dex1, int dex2); /* square matrices only */

	/* perform Gaussian elimination with the given RHS vector */
	/* Note: matrix overwritten during solution!              */
	void LinearSolve(dArrayT& RHS);
	void LinearSolve2(dArrayT& RHS); /* skips zeroes */

	/* assignment operators */
	LAdMatrixT& operator=(const dMatrixT& RHS);
	LAdMatrixT& operator=(const double value);
};

/* inlines */

/* assignment operators */
inline LAdMatrixT& LAdMatrixT::operator=(const dMatrixT& RHS)
{
	/* inherited */
	dMatrixT::operator=(RHS);
	return *this;
}

inline LAdMatrixT& LAdMatrixT::operator=(const double value)
{
	/* inherited */
	dMatrixT::operator=(value);
	return *this;
}

#endif /* _LA_DMATRIX_T_H_ */
