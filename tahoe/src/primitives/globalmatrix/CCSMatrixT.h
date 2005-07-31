/* $Id: CCSMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (05/29/1996)                                          */
/* This is the interface for a Symmetric matrix stored in                 */
/* Compact Column form.                                                   */
/* To initialize:                                                         */
/* 			(1) call constructor with system dimension dim                */
/* 			(2) set the ColumnHeight for column = 0...dim-1               */
/* 			(3) call Initialize() to allocate space for the               */
/* matrix and set the diagonal position array.                            */
/* Note: assembly positions (equation numbers) = 1...fNumEQ               */

#ifndef _CCSMATRIX_T_H_
#define _CCSMATRIX_T_H_

/* project headers */
#include "Environment.h"
#include "ExceptionCodes.h"

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "LinkedListT.h"

class CCSMatrixT: public GlobalMatrixT
{
public:

	/* constructor */
	CCSMatrixT(ostream& out, int check_code);

	/* destructor */	
	virtual ~CCSMatrixT(void);
	
	/* set the internal matrix structure.
	 * NOTE: do not call Initialize() equation topology has been set
	 * with AddEquationSet() for all equation sets */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* set all matrix values to 0.0 */
	virtual void Clear(void);
		
	/* add element group equations to the overall topology.
	 * NOTE: assembly positions (equation numbers) = 1...fNumEQ
	 * equations can be of fixed size (iArray2DT) or
	 * variable length (RaggedArray2DT) */
	virtual void AddEquationSet(const iArray2DT& eqset);
	virtual void AddEquationSet(const RaggedArray2DT<int>& eqset);
	
	/* assemble the element contribution into the LHS matrix - assumes
	 * that elMat is square (n x n) and that eqnos is also length n.
	 * NOTE: assembly positions (equation numbers) = 1...fNumEQ */
	virtual void Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos);
	
	/* compute the diagonal weighted residual norm - no check as
	 * to whether the matrix is factorized or not */
	double ResidualNorm(const dArrayT& result) const;
	
	/* compute the sum of the elements on the prescribed row/col,
	 * where rownum = 0...fNumEQ-1 */
	double AbsRowSum(int rownum) const;
	
	/* multiply the matrix with d and return the result in Kd.
	 * Note: do not call if the matrix is already factorized */
	void MultKd(const dArrayT& d, dArrayT& Kd) const;

	/* return the value of p_i K_ij p_j */
	double pTKp(const dArrayT& p) const;

	/* returns 1 if the factorized matrix contains a negative
	 * pivot.  Matrix MUST be factorized.  Otherwise function
	 * returns 0 */
	int HasNegativePivot(void) const;

	/* assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);

	/* TESTING: write non-zero elements of matrix in Aztec readable
	 *          format */
	void WriteAztecFormat(ostream& out) const;

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

protected:

	/* element accessor - READ ONLY */
	double operator()(int row, int col) const;			
	
	/* output operator */
	friend ostream& operator<<(ostream& out, const CCSMatrixT& matrix);	

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

	/* rank check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;

	/* Returns the number of elements ABOVE the diagonal in col */
	int ColumnHeight(int col) const;

private:

	/* (re-) compute the matrix structure and return the bandwidth
	 * and mean bandwidth */
	void ComputeSize(int& num_nonzero, int& mean_bandwidth, int& bandwidth);

	/* computes the column heights for the given equation list */
	void SetColumnHeights(const iArray2DT& eqnos);
	void SetColumnHeights(const RaggedArray2DT<int>& eqnos);

	/* Returns number of non-zero elements */
	int  NumberOfFilled(int printRCV = 0);
	void FillWithOnes(const iArray2DT& eqnos);
	void FillWithOnes(const RaggedArray2DT<int>& eqnos);
		
protected:

	/* equations sets */
	LinkedListT<const iArray2DT*> fEqnos;
	LinkedListT<const RaggedArray2DT<int>*> fRaggedEqnos;

	int*    fDiags; 	
	int     fNumberOfTerms;
	double* fMatrix;
};

/* Returns the number of elements ABOVE the diagonal in col */
inline int CCSMatrixT::ColumnHeight(int col) const
{
#if __option (extended_errorcheck)
	if (col < 0 || col >= fLocNumEQ) throw eGeneralFail;
#endif

	return ( (col > 0) ? (fDiags[col] - fDiags[col-1] - 1) : 0);
}

#endif /* _CCSMATRIX_T_H_ */
