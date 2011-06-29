/* $Id: AztecMatrixT.h,v 1.16 2005-04-13 21:50:27 paklein Exp $ */
/* created: paklein (08/10/1998) */
#ifndef _AZTEC_MATRIX_T_H_
#define _AZTEC_MATRIX_T_H_

/* base classes */
#include "GlobalMatrixT.h"

/* library support option */
#ifdef __AZTEC__

/* direct members */
#include "AutoArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class Aztec_fe;
class ifstreamT;
class ParameterListT;

/** interface for Aztec linear solver library */
class AztecMatrixT: public GlobalMatrixT
{
public:

	/* constuctor */
	AztecMatrixT(ostream& out, int check_code, const CommunicatorT& comm,
		const ParameterListT& parameters);

	/** copy constructor */
	AztecMatrixT(const AztecMatrixT& source);
	
	/* destuctor */
	virtual ~AztecMatrixT(void);

	/* set matrix structure and allocate space.
	 * NOTE: do not call Initialize until all equations sets have been
	 * registered with SetStructure */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** write information to output stream after AztecMatrixT::Initialize
	 * has been called */
	virtual void Info(ostream& out);
	
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
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operator */
	AztecMatrixT& operator=(const AztecMatrixT& rhs);

	/** return a clone of self. Caller is responsible for disposing of the matrix */
	virtual GlobalMatrixT* Clone(void) const;

protected:
	
	/* determine new search direction and put the results in result */
	virtual void BackSubstitute(dArrayT& result);

	/* rank check functions */
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(bool force = false) const;

private:

	/* C++ interface to Aztec 1.0 */
	Aztec_fe* fAztec;
	
	/* for assembly operations */
	AutoArrayT<int> fRowDexVec, fColDexVec;
	AutoArrayT<double> fValVec;
	dMatrixT fValMat;
};

/* library support options */
} // namespace Tahoe 
#endif /* __AZTEC__ */
#endif /* _AZTEC_MATRIX_T_H_ */
