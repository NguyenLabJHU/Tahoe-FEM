/* $Id: CCSMatrixT.h,v 1.15 2004-10-04 18:40:51 paklein Exp $ */
/* created: paklein (05/29/1996) */
#ifndef _CCSMATRIX_T_H_
#define _CCSMATRIX_T_H_

/* project headers */
#include "Environment.h"
#include "ExceptionT.h"

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "LinkedListT.h"

namespace Tahoe {

/** This is the interface for a Symmetric matrix stored in
 * Compact Column form.
 * To initialize:
 * -# call constructor with system dimension dim
 * -# set the ColumnHeight for column = 0...dim-1
 * -# call Initialize() to allocate space for the matrix and set the diagonal position array.
 * Note: assembly positions (equation numbers) = 1...fNumEQ
 */
class CCSMatrixT: public GlobalMatrixT
{
public:

	/** constructor */
	CCSMatrixT(ostream& out, int check_code);

	/** copy constructor */
	CCSMatrixT(const CCSMatrixT& source);

	/** destructor */	
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
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);
	
	/* compute the diagonal weighted residual norm - no check as
	 * to whether the matrix is factorized or not */
	double ResidualNorm(const dArrayT& result) const;
	
	/* compute the sum of the elements on the prescribed row/col,
	 * where rownum = 0...fNumEQ-1 */
	double AbsRowSum(int rownum) const;

	/* return the value of p_i K_ij p_j */
	double pTKp(const dArrayT& p) const;

	/* returns 1 if the factorized matrix contains a negative
	 * pivot.  Matrix MUST be factorized.  Otherwise function
	 * returns 0 */
	int HasNegativePivot(void) const;

	/* TESTING: write non-zero elements of matrix in Aztec readable
	 *          format */
	void WriteAztecFormat(ostream& out) const;

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kSymmetric; };

	/* find the smallest and largest diagonal value */
	void FindMinMaxPivot(double& min, double& max, double& abs_min, double& abs_max) const;

	/** assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& rhs);
	
	/** return a clone of self. Caller is responsible for disposing of the matrix */
	virtual GlobalMatrixT* Clone(void) const;

	/** matrix-vector product. Works only if called before the matrix has been
	 * factorized. 
	 * \param x vector to use for calculating the product
	 * \param b destination for the result 
	 * \return true if the product was calculated successful */
	virtual bool Multx(const dArrayT& x, dArrayT& b) const;

	/** Tranpose[matrix]-vector product. Works only if called before the matrix 
	 * has been factorized. 
	 * \param x vector to use for calculating the product
	 * \param b destination for the result
	 * \return true if the product was calculated successful */
	virtual bool MultTx(const dArrayT& x, dArrayT& b) const;

	/** return the values along the diagonal of the matrix. Derived classes
	 * must reimplement this function to extrat the diagonals from the
	 * matrix-specific storage schemes.
	 * \param diags returns with the diagonals of the matrix if the function
	 *        is supported. Otherwise is left unchanged.
	 * \return true if the diagonal values where collected successfully */
	virtual bool CopyDiagonal(dArrayT& diags) const;

	/** \name check functions */
	/*@{*/
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(bool force = false) const;
	/*@}*/

protected:

	/** \name element accessors */
	/*@{*/
	/** read only accessor returns 0.0 for nonstored values */
	double operator()(int row, int col) const;

	/** read/write accessor throws exception for nonstored values */
	double& operator()(int row, int col);
	/*@}*/

	/* output operator */
	friend ostream& operator<<(ostream& out, const CCSMatrixT& matrix);	

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

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

	/** \name equations sets */
	/*@{*/
	LinkedListT<const iArray2DT*> fEqnos;
	LinkedListT<const RaggedArray2DT<int>*> fRaggedEqnos;
	/*@}*/

	/** \name the matrix */
	/*@{*/
	int*    fDiags; 	
	int     fNumberOfTerms;
	double* fMatrix;
	/*@}*/
	
	/** runtime flag */
	bool fIsFactorized;
};

/* Returns the number of elements ABOVE the diagonal in col */
inline int CCSMatrixT::ColumnHeight(int col) const
{
#if __option (extended_errorcheck)
	if (col < 0 || col >= fLocNumEQ) throw ExceptionT::kGeneralFail;
#endif

	return ( (col > 0) ? (fDiags[col] - fDiags[col-1] - 1) : 0);
}

/* Tranpose[matrix]-vector product */
inline bool CCSMatrixT::MultTx(const dArrayT& x, dArrayT& b) const
{
	return CCSMatrixT::Multx(x, b);
}

/* element accessor */
inline double& CCSMatrixT::operator()(int row, int col)
{
	const char caller[] = "CCSMatrixT::operator()";

#if __option(extended_errorcheck)
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) ExceptionT::OutOfRange(caller);
	if (col < 0 || col >= fLocNumEQ) ExceptionT::OutOfRange(caller);
#endif
	
	if (row == col) /* element on the diagonal */
		return fMatrix[fDiags[col]];
	else
	{
		/* look into upper triangle */
		int& r = (row > col) ? col : row;
		int& c = (row > col) ? row : col;

		int colht = ColumnHeight(c);
		int hrow  = c - r;		
		if (hrow > colht) /* element above the skyline */
			ExceptionT::OutOfRange(caller);
		
		return fMatrix[fDiags[c] - hrow];
	}
}

} /* namespace Tahoe */

#endif /* _CCSMATRIX_T_H_ */
