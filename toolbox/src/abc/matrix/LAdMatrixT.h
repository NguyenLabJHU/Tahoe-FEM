/* $Id: LAdMatrixT.h,v 1.5 2003-11-21 22:41:36 paklein Exp $ */
/* created: paklein (12/05/1996)                                          */
/* dMatrixT with some linear algebra functions                            */

#ifndef _LA_DMATRIX_T_H_
#define _LA_DMATRIX_T_H_

/* base class */
#include "dMatrixT.h"

namespace Tahoe {

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

//*****************************************************************
// Iterative  routine -- BiCGSTAB
//
// BiCGSTAB solves the unsymmetric linear system Ax = RHS 
// using the Preconditioned BiConjugate Gradient Stabilized method
//
// BiCGSTAB follows the algorithm described on p. 27 of the 
// SIAM Templates book.
//
//
// Input Values
//        x  -- initial guess to Ax = RHS
//      RHS  -- the right hand side of the system Ax = RHS
//		  M  -- Diagonal preconditioner
//				double 	M ->   M*I,  ( I is the identity matrix )
//				dArrayT M ->   M,    ( M is a a digonal matrix with M(i) -> M(i,i) )
//
// 	 max_it	 --  Maximum iterations before failure. Should be on the order of Dim(A).  
//      tol  --  the residual below which the solution is considered to be convereged.
//				 tol = |residual|/|RHS|
//               Should be on the order of machine precision
//   
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//	return values:
//		   0 -- convergence within max_iter
//		   1 -- did not converge to the specified tolerance within max_iter
//		   2 -- internal breakdown
//		   3 -- internal breakdown
//  
//*****************************************************************
	int BiCGStab(dArrayT& x ,const dArrayT& RHS ,const double   M , int max_it , double tol );
	int BiCGStab(dArrayT& x ,const dArrayT& RHS ,const dArrayT& M , int max_it , double tol );

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

} // namespace Tahoe 
#endif /* _LA_DMATRIX_T_H_ */
