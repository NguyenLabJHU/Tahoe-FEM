/* $Id: FullMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (03/07/1998)                                          */
/* Virtual base class for all global matrix objects                       */

#ifndef _FULL_MATRIX_T_H_
#define _FULL_MATRIX_T_H_

/* base class */
#include "GlobalMatrixT.h"

/* direct members */
#include "LAdMatrixT.h"

class FullMatrixT: public GlobalMatrixT
{
public:

	/* constructor */
	FullMatrixT(ostream& out, int check_code);
		
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

	/* strong manipulation functions */
	//TEMP should be pure virtual, but no time to update others
	//     so just throw exception for now
	virtual void OverWrite(const ElementMatrixT& elMat, const iArrayT& eqnos);
	virtual void Disassemble(dMatrixT& elMat, const iArrayT& eqnos) const;

	/* assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;
	
protected:

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* solution driver */
	virtual void BackSubstitute(dArrayT& result);

	/* rank check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;
	
protected:

	/* the matrix */
	LAdMatrixT fMatrix;
};

#endif /* _FULL_MATRIX_T_H_ */
