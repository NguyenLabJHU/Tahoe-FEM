/* $Id: DiagonalMatrixT.h,v 1.16 2005-02-04 22:01:54 paklein Exp $ */
/* created: paklein (03/23/1997) */

#ifndef _DIAGONAL_MATRIX_H_
#define _DIAGONAL_MATRIX_H_

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "dArrayT.h"

namespace Tahoe {

/** diagonal matrix */
class DiagonalMatrixT: public GlobalMatrixT
{
public:

	/** enum to signal how to assemble non-diagonal contributions to the matrix */
	enum AssemblyModeT {kNoAssembly = 0, /**< do not assemble, throw ExceptionT::xception */ 
	                    kDiagOnly   = 1, /**< assemble the diagonal values only */
                        kAbsRowSum  = 2  /**< assemble the L1 norm of the row */};

	/** constructors */
	DiagonalMatrixT(ostream& out, int check_code, AssemblyModeT mode);

	/** copy constructor */
	DiagonalMatrixT(const DiagonalMatrixT& source);

	/** DiagonalMatrixT::Solve does preserve the data in the matrix */
	virtual bool SolvePreservesData(void) const { return true; };	  

	/* set assemble mode */
	void SetAssemblyMode(AssemblyModeT mode);

	/* set the internal matrix structure.
	 * NOTE: do not call Initialize() equation topology has been set
	 * with AddEquationSet() for all equation sets */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* set all matrix values to 0.0 */
	virtual void Clear(void);
		
	/* add element group equations to the overall topology.
	 * NOTE: assembly positions (equation numbers) = 1...fDimension
	 * equations can be of fixed size (iArray2DT) or
	 * variable length (RaggedArray2DT) */
	virtual void AddEquationSet(const iArray2DT& eqset);
	virtual void AddEquationSet(const RaggedArray2DT<int>& eqset);
	
	/* assemble the element contribution into the LHS matrix - assumes
	 * that elMat is square (n x n) and that eqnos is also length n.
	 * NOTE: assembly positions (equation numbers) = 1...fDimension */
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);

	/* fetch values */
	virtual void DisassembleDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const;

	/* access to the data */
	dArrayT& TheMatrix(void);

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kDiagonal; };

	/** assignment operator */
	DiagonalMatrixT& operator=(const DiagonalMatrixT& rhs);

	/** return a clone of self. Caller is responsible for disposing of the matrix */
	virtual GlobalMatrixT* Clone(void) const;

	/** matrix-vector product. OK to call either before or after the matrix is
	 * factorized */
	virtual void Multx(const dArrayT& x, dArrayT& b) const;

	/** Tranpose[matrix]-vector product. OK to call either before or after the matrix 
	 * is factorized */
	virtual void MultTx(const dArrayT& x, dArrayT& b) const;

	/** vector-matrix-vector product */
	virtual double MultmBn(const dArrayT& m, const dArrayT& n) const;
	
protected:

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* solution driver */
	virtual void BackSubstitute(dArrayT& result);

	/* check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(bool force = false) const;
	
private:

	/** the matrix */
	dArrayT	fMatrix;

	/** runtime flag */
	bool fIsFactorized;

	/* mode flag */
	AssemblyModeT fMode;
};

/* inlines */

/* access to the data */
inline dArrayT& DiagonalMatrixT::TheMatrix(void)
{
	return fMatrix;
}

/* Tranpose[matrix]-vector product */
inline void DiagonalMatrixT::MultTx(const dArrayT& x, dArrayT& b) const { 
	DiagonalMatrixT::Multx(x, b);
}; 

} // namespace Tahoe 
#endif /* _DIAGONAL_MATRIX_H_ */
