/* $Id: GlobalMatrixT.h,v 1.16 2005-01-07 21:22:49 paklein Exp $ */
/* created: paklein (03/23/1997) */
#ifndef _GLOBAL_MATRIX_H_
#define _GLOBAL_MATRIX_H_

#include "GlobalT.h"
#include "Environment.h"
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class dMatrixT;
class ElementMatrixT;
template <class TYPE> class ArrayT;
template <class nTYPE> class nArrayT;
class dArrayT;
class iArray2DT;
template <class TYPE> class RaggedArray2DT;

/** virtual base class for global matrix objects */
class GlobalMatrixT
{
public:

	/** check codes */
	enum CheckCodeT {kNoCheck = 0,
                  kZeroPivots = 1,
                   kAllPivots = 2,
                    kPrintLHS = 3,
                    kPrintRHS = 4,
               kPrintSolution = 5,
                    kCheckLHS = 6};

	/** equation numbering scope */
	enum EquationNumberScopeT {
		kLocal  = 0,
		kGlobal = 1}; // for parallel solvers

	/** constructor */
	GlobalMatrixT(ostream& out, int check_code);

	/** copy constructor */
	GlobalMatrixT(const GlobalMatrixT& source);

	/** destructor */	
	virtual ~GlobalMatrixT(void);

	/** return true if GlobalMatrixT::Solve preserves the data in the matrix or
	 * false if the solution procedure overwrites this data. Preserving the data
	 * implies that all of the operations in the matrix are valid both before and
	 * after a call to GlobalMatrixT::Solve even if the matrix is not refilled. */
	virtual bool SolvePreservesData(void) const { return false; };	  

	/** set the internal matrix structure.
	 * NOTE: do not call Initialize() equation topology has been set
	 * with AddEquationSet() for all equation sets */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/** clear values for next assembly */
	virtual void Clear(void) {};
	
	/** solve for rhs passed in result and overwritten with solution */
	bool Solve(dArrayT& result);
	
	/** \name add element group equations to the overall topology
	 * Assembly positions (equation numbers) = 1...fDimension
	 * equations can be of fixed size (iArray2DT) or
	 * variable length (RaggedArray2DT) */
	/*@{*/
	virtual void AddEquationSet(const iArray2DT& eqset) = 0;
	virtual void AddEquationSet(const RaggedArray2DT<int>& eqset) = 0;
	/*@}*/

	/** \name assemble operators
	 * Assemble the element contribution into the LHS matrix. Assumes
	 * that elMat is square (n x n) and that eqnos is also length n.
	 * NOTE: assembly positions (equation numbers) = 1...fDimension */
	/*@{*/
	/** assembly of square element matrix. The global equation numbers associated
	 * with the rows and columns of the matrix are the same. */
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos) = 0;

	/** assembly of general element matrix. The global equation numbers associated
	 * with the rows and columns of the matrix are specified separately and the
	 * matrix does not need to be square. */
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos) = 0;

	/** assembly of a diagonal matrix */
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos) = 0;
	/*@}*/

	/* strong manipulation functions 
	 * NOTE: These must be overridden to provide support for these functions.
	 *       By default, these all throw ExceptionT::xceptions. These could be pure
	 *       virtual, but that requires updating all derived matrix types */
	/*@{*/
	virtual void OverWrite(const ElementMatrixT& elMat, const nArrayT<int>& eqnos);
	virtual void Disassemble(dMatrixT& matrix, const nArrayT<int>& eqnos) const;
	virtual void DisassembleDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const;
	/*@}*/

	/** \name number scope and reordering */
	/*@{*/
	virtual EquationNumberScopeT EquationNumberScope(void) const = 0;
	virtual bool RenumberEquations(void) const = 0;
	/*@}*/
	
	/** \name accessors */
	/*@{*/
	int CheckCode(void) const;
	int NumEquations(void) const;
	int StartEquation(void) const;
	
	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const = 0;
	/*@}*/

	/** assignment operator */
	GlobalMatrixT& operator=(const GlobalMatrixT& rhs);
	
	/** return a clone of self. Caller is responsible for disposing of the matrix */
	virtual GlobalMatrixT* Clone(void) const = 0;

	/** matrix-vector product. Derived classes should reimplement this
	 * function if the product is supported. 
	 * \param x vector to use for calculating the product
	 * \param b destination for the result 
	 * \return true if the product was calculated successful */
	virtual bool Multx(const dArrayT& x, dArrayT& b) const;

	/** Tranpose[matrix]-vector product. Derived classes should reimplement this
	 * function if the product is supported.
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
	virtual void PrintAllPivots(void) const = 0;
	virtual void PrintZeroPivots(void) const = 0;

	/** write matrix if check code is GlobalMatrixT::kPrintLHS or if force is true */
	virtual void PrintLHS(bool force = false) const = 0;
	/*@}*/	

protected:

	/** precondition matrix */
	virtual void Factorize(void) {};
	
	/** solution driver */
	virtual void BackSubstitute(dArrayT& result) = 0;

	/** \name check functions */
	/*@{*/
	void PrintRHS(const dArrayT& RHS) const;
	void PrintSolution(const dArrayT& solution) const;
	/*@}*/

	/** \name inline functions needed during factorization */
	/*@{*/
	static int Min(int a, int b);
	static int Max(int a, int b);
	static double Max(double a, double b);
	static double Min(double a, double b);
	static double Dot(double* vec1, double* vec2, int length);	
	/*@}*/

protected:

	/** output stream */
	ostream& fOut;

	/** \name parameters */
	/*@{*/
	int fCheckCode;  	
	int	fTotNumEQ;
	int	fLocNumEQ;
	int fStartEQ; //1,...
	/*@}*/
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

} // namespace Tahoe 
#endif /* _GLOBAL_MATRIX_H_ */
