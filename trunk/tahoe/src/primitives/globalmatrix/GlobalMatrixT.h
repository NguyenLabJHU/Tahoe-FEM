/* $Id: GlobalMatrixT.h,v 1.1.1.1 2001-01-29 08:20:23 paklein Exp $ */
/* created: paklein (03/23/1997)                                          */
/* Virtual base class for all global matrix objects                       */

#ifndef _GLOBAL_MATRIX_H_
#define _GLOBAL_MATRIX_H_

#include "Environment.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class dMatrixT;
class ElementMatrixT;
class iArrayT;
class dArrayT;
class iArray2DT;
template <class TYPE> class RaggedArray2DT;

class GlobalMatrixT
{
public:

	/* check codes */
	enum CheckCodeT {kNoCheck = 0,
                  kZeroPivots = 1,
                   kAllPivots = 2,
                    kPrintLHS = 3,
                    kPrintRHS = 4,
               kPrintSolution = 5};

	/* equation numbering scope */
	enum EquationNumberScopeT {
		kLocal  = 0,
		kGlobal = 1}; // for parallel solvers

	/* constructor */
	GlobalMatrixT(ostream& out, int check_code);

	/* destructor */	
	virtual ~GlobalMatrixT(void);

	/* set the internal matrix structure.
	 * NOTE: do not call Initialize() equation topology has been set
	 * with AddEquationSet() for all equation sets */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* set all matrix values to 0.0 */
	virtual void Clear(void) = 0;
	
	/* solve for rhs passed in result and overwritten with solution */
	void Solve(dArrayT& result);
	
	/* add element group equations to the overall topology.
	 * NOTE: assembly positions (equation numbers) = 1...fDimension
	 * equations can be of fixed size (iArray2DT) or
	 * variable length (RaggedArray2DT) */
	virtual void AddEquationSet(const iArray2DT& eqset) = 0;
	virtual void AddEquationSet(const RaggedArray2DT<int>& eqset) = 0;
	
	/* assemble the element contribution into the LHS matrix - assumes
	 * that elMat is square (n x n) and that eqnos is also length n.
	 * NOTE: assembly positions (equation numbers) = 1...fDimension */
	virtual void Assemble(const ElementMatrixT& elMat, const iArrayT& eqnos) = 0;

	/* strong manipulation functions */
	virtual void OverWrite(const ElementMatrixT& elMat, const iArrayT& eqnos);
	virtual void Disassemble(dMatrixT& elMat, const iArrayT& eqnos) const;
		//TEMP should be pure virtual, but no time to update others
		//     so just throw exception for now

	/* assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& RHS);
	
	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const = 0;
	virtual bool RenumberEquations(void) const = 0;
	
	/* accessors */
	int CheckCode(void) const;
	int NumEquations(void) const;
	int StartEquation(void) const;
	
protected:

	/* precondition matrix */
	virtual void Factorize(void) = 0;
	
	/* solution driver */
	virtual void BackSubstitute(dArrayT& result) = 0;

	/* check functions */
	virtual void PrintAllPivots(void) const = 0;
	virtual void PrintZeroPivots(void) const = 0;
	virtual void PrintLHS(void) const = 0;
	void PrintRHS(const dArrayT& RHS) const;
	void PrintSolution(const dArrayT& solution) const;

	/* inline functions needed during factorization */
	static int Min(int a, int b);
	static int Max(int a, int b);
	static double Max(double a, double b);
	static double Min(double a, double b);
	static double Dot(double* vec1, double* vec2, int length);	
	
protected:

	/* output stream */
	ostream& fOut;

	/* parameters */
	int fCheckCode;  	
	int	fTotNumEQ;
	int	fLocNumEQ;
	int fStartEQ; //1,...
	
	/* runtime flag */
	int fIsFactorized;
};

/* return the check code */
inline int GlobalMatrixT::CheckCode(void) const { return fCheckCode; }
inline int GlobalMatrixT::NumEquations(void) const { return fLocNumEQ; }
inline int GlobalMatrixT::StartEquation(void) const { return fStartEQ; }

/* Inline functions: Protected */
inline int GlobalMatrixT::Min(int a, int b)
{	
	return (a > b) ? b : a;
}

inline int GlobalMatrixT::Max(int a, int b)
{	
	return (a < b) ? b : a;
}

inline double GlobalMatrixT::Min(double a, double b)
{	
	return (a > b) ? b : a;
}

inline double GlobalMatrixT::Max(double a, double b)
{	
	return (a < b) ? b : a;
}

inline double GlobalMatrixT::Dot(double* vec1, double* vec2, int length)
{
	register double dot = 0.0;	
	for (int i = 0; i < length; i++)
		dot += (*vec1++)*(*vec2++);
	return dot;
}

#endif /* _GLOBAL_MATRIX_H_ */
