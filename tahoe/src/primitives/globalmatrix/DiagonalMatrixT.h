/* $Id: DiagonalMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Virtual base class for all global matrix objects                       */

#ifndef _DIAGONAL_MATRIX_H_
#define _DIAGONAL_MATRIX_H_

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "dArrayT.h"

class DiagonalMatrixT: public GlobalMatrixT
{
public:

	enum AssemblyModeT {kNoAssembly = 0,
	                    kDiagOnly   = 1,
                        kAbsRowSum  = 2};

	/* constructors */
	DiagonalMatrixT(ostream& out, int check_code, AssemblyModeT mode);

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
	virtual void Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos);

	/* assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);

	/* access to the data */
	dArrayT& TheMatrix(void);

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;
	
protected:

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* solution driver */
	virtual void BackSubstitute(dArrayT& result);

	/* check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;
	
private:

	dArrayT	fMatrix;
	
	/* mode flag */
	AssemblyModeT fMode;
};

/* inlines */

/* access to the data */
inline dArrayT& DiagonalMatrixT::TheMatrix(void)
{
	return fMatrix;
}

#endif /* _DIAGONAL_MATRIX_H_ */
