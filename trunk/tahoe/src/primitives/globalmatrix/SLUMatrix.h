/* $Id: SLUMatrix.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: rbridson (06/30/2000)                                         */
/* Interface to SuperLU solver library, extending GlobalMatrixT           */

#ifndef _SLU_MATRIX_H_
#define _SLU_MATRIX_H_

/* project headers */
#include "Environment.h"

/* library support */
#ifdef __SUPERLU__

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "LinkedListT.h"

/* external SuperLU stuff */
#include "superlu.h"

class SLUMatrix: public GlobalMatrixT
{
public:

	/* constructor */
	SLUMatrix(ostream& out, int check_code);

	/* destructor */	
	~SLUMatrix(void);

	/* set the internal matrix structure.
	 * NOTE: do not call Initialize() until equation topology has been set
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

	/* assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);
	
	/* element accessor - READ ONLY */
	double Element(int row, int col) const;

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

protected:

	/* constants to be tuned for efficiency of factorization */
	static const int kSupernodeRelax = 5;
	static const int kPanelSize = 10;

	/* output in sparse format */
	friend ostream& operator<<(ostream& out, const SLUMatrix& matrix);

	/* decompose matrix into PLU */
	virtual void Factorize(void);
	
	/* solution driver */
	virtual void BackSubstitute(dArrayT& result);

	/* check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;

	/* element accessor - read and write, for assembly. Exception for */
	/* access to unstored zero. */
	double& operator()(int row, int col) const;

	/* (over)estimate number of nonzeros from equation sets */
	virtual void EstimateNNZ (int *colLength, int &nnz);

	/* insert all the element equations into A */
	virtual void InsertEquations (NCformat *A, int *colLength, int &nnz);

	/* Insert the list nzlist of nonzeros (with length nzlen) into column c
	 * of matrix A, using colLengths to keep track of column lengths, and
	 * keeping nnz up-to-date. The columns of A will have nonzeros in
	 * ascending order. */
	virtual void InsertNZ (NCformat *A, int *colLength, int &nnz, int c,
	   int nzlen, int *nzlist);

	/* compress columns in A */
	virtual void CompressColumns (NCformat *A, const int *colLength);

	/* figure out a good ordering of the columns of the matrix in perm_c */
	virtual void OrderColumns (void);

protected:

	/* The matrix and its factors in SuperLU column format */
	SuperMatrix fMatrix;
	SuperMatrix fLower;
	SuperMatrix fUpper;

	/* flag to see if fLower and fUpper have been allocated */
	bool fLUallocated;

	/* Column and row permutations used in factorization */
	int *fPerm_c, *fPerm_r;
	bool fIsColOrdered;   // true if perm_c has been set

	/* Symbolic information used in factorization */
	int *fEtree;

	/* equation sets */
	LinkedListT<const iArray2DT*> fEqnos;
	LinkedListT<const RaggedArray2DT<int>*> fRaggedEqnos;
};

#endif /* __SUPERLU__ */
#endif /* _SLU_MATRIX_H_ */
