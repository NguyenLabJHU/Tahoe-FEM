/* $Id: AztecMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (08/10/1998)                                          */
/* interface for using Aztec with fe++                                    */

#ifndef _AZTEC_MATRIX_T_H_
#define _AZTEC_MATRIX_T_H_

/* base classes */
#include "GlobalMatrixT.h"

/* library support option */
#ifdef __AZTEC__

/* direct members */
#include "AutoArrayT.h"
#include "dMatrixT.h"

/* foward declarations */
class Aztec_fe;
class ifstreamT;

class AztecMatrixT: public GlobalMatrixT
{
public:

	/* constuctor */
	AztecMatrixT(ifstreamT& in, ostream& out, int check_code);

	/* destuctor */
	virtual ~AztecMatrixT(void);

	/* set matrix structure and allocate space.
	 * NOTE: do not call Initialize until all equations sets have been
	 * registered with SetStructure */
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
	 *
	 * NOTE: assembly positions (equation numbers) = 1...fNumEQ */
	virtual void Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos);

	/* assignment operator - not implemented */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

protected:

	/* precondition matrix */
	virtual void Factorize(void);
	
	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

	/* rank check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(void) const;

private:

	/* input stream */
	ifstreamT& fInput;
	
	/* C++ interface to Aztec 1.0 */
	Aztec_fe* fAztec;
	
	/* for assembly operations */
	AutoArrayT<int> fRowDexVec, fColDexVec;
	AutoArrayT<double> fValVec;
	dMatrixT fValMat;
};

/* library support options */
#endif /* __AZTEC__ */
#endif /* _AZTEC_MATRIX_T_H_ */
