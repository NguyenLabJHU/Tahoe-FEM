/* $Id: SuperLU_DISTMatrixT.h,v 1.1 2004-03-16 10:03:21 paklein Exp $ */
#ifndef _SUPER_LU_DIST_MATRIX_T_H_
#define _SUPER_LU_DIST_MATRIX_T_H_

//TEMP
//#define __SUPERLU_DIST__
//#define __TAHOE_MPI__

/* library support */
#ifdef __SUPERLU_DIST__
#ifdef __TAHOE_MPI__

/* base class */
#include "GlobalMatrixT.h"

/* SuperLU type definitions */
#include "superlu_ddefs.h"

namespace Tahoe {

/* forward declarations */
class MSRBuilderT;
class CommunicatorT;

/** interface to SuperLU 3.0 serial linear solver */
class SuperLU_DISTMatrixT: public GlobalMatrixT
{
public:

	/** constructor */
	SuperLU_DISTMatrixT(ostream& out, int check_code, bool symmetric, CommunicatorT& comm);

	/** destructor */
	~SuperLU_DISTMatrixT(void);

	/** set the internal matrix structure */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* set all matrix values to 0.0 */
	virtual void Clear(void);
	
	/** \name set matrix structure
	 * active equations are eq > 0 */
	/*@{*/
	void AddEquationSet(const iArray2DT& eqnos);
	void AddEquationSet(const RaggedArray2DT<int>& eqnos);
	/*@}*/

	/** \name assembly operators
	 * active equations are eq > 0 */	
	/*@{*/
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& eqnos);
	virtual void Assemble(const ElementMatrixT& elMat, const ArrayT<int>& row_eqnos,
		const ArrayT<int>& col_eqnos);
	virtual void Assemble(const nArrayT<double>& diagonal_elMat, const ArrayT<int>& eqnos);
	/*@}*/
	
	/* element accessor - READ ONLY */
	double Element(int row, int col) const;

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operator */
	virtual GlobalMatrixT& operator=(const GlobalMatrixT& rhs);

	/** return a clone of self */
	virtual GlobalMatrixT* Clone(void) const;

protected:

	/** \name tuning parameters */
	/*@{*/
	static const int kSupernodeRelax = 5;
	static const int kPanelSize = 10;
	/*@}*/

	/** output in sparse format */
	friend ostream& operator<<(ostream& out, const SuperLU_DISTMatrixT& matrix) {
		NCformat *A = (NCformat*) (matrix.fA.Store);
		int i, j;
		for (i = 0; i < matrix.fA.ncol; ++i) {
			out << -1 << " " << 0.0 << "\n";
			for (j = A->colptr[i]; j < A->colptr[i+1]; ++j) {
				out << A->rowind[j]+1 << " ";
				out << ((double*)A->nzval)[j] << "\n";
			}
		}
		return out;
	};

	/** solution driver. Calls all-in-one driver provided with SuperLU 3.0 which 
	 * can be called for solving multiple right-hand sides or just resolving
	 * a matrix with the same sparsity pattern as a previous solve. This driver
	 * routine is adapted from dlinsolx2.c provided in the SuperLU 3.0 examples. */
	virtual void BackSubstitute(dArrayT& result);

	/** \name check functions */
	/*@{*/
	virtual void PrintAllPivots(void) const;
	virtual void PrintZeroPivots(void) const;
	virtual void PrintLHS(bool force) const;
	/*@}*/

	/** element accessor - read and write, for assembly. Exception for */
	/* access to unstored zero. */
	double& operator()(int row, int col);

private:

	/** no copy constructor */
	SuperLU_DISTMatrixT(const SuperLU_DISTMatrixT&);

protected:

	/** MP environment */
	CommunicatorT& fComm;	

	/** matrix structure builder */
	MSRBuilderT* fBuilder;

	/** \name matrix and factors in SuperLU formats */
	/*@{*/
	SuperMatrix fA;
	AutoArrayT<int> frowptr;
	AutoArrayT<int> fcolind;
	AutoArrayT<double> fnzval;

	SuperMatrix fL;
	SuperMatrix fU;
	/*@}*/

	/** \name factorization workspace */
	/*@{*/
	/** column permutations */
	AutoArrayT<int> fperm_c;

	/** row permutations */
	AutoArrayT<int> fperm_r;

	/** symbolic information used in factorization */
	AutoArrayT<int> fetree;

	/** rhs vector used for linear solves */
	SuperMatrix fB;

	/** solution vector */
	SuperMatrix fX;
	
	/** row scaling */
	AutoArrayT<double> fR;

	/** column scaling */
	AutoArrayT<double> fC;
	
	/** equilibration */
	char fequed;
	/*@}*/

	/** SuperLU options */
	superlu_options_t foptions;

	/** \name factorization flags */
	/*@{*/
	/** true if matrix has been symbolically factorized */
	bool fIsSymFactorized;

	/** true if matrix has been numerically factorized */
	bool fIsNumFactorized;
	/*@}*/
};

/* element accessor - read and write, for assembly. Exception for */
/* access to unstored zero. row and col are */
inline double& SuperLU_DISTMatrixT::operator()(int row, int col)
{
	const char caller[] = "SuperLU_DISTMatrixT::operator()";

#if __option(extended_errorcheck)
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) ExceptionT::GeneralFail(caller);
	if (col < 0 || col >= fLocNumEQ) ExceptionT::GeneralFail(caller);
#endif

	NCformat *A = (NCformat*) fA.Store;

	/* look through column col for the given row index */
	int start = A->colptr[col];
	int stop  = A->colptr[col+1];
	for (int i = start; i < stop; ++i)
	{
		if (row == A->rowind[i])
			return ((double*)A->nzval)[i];
	}

	/* otherwise this nonzero wasn't present */
	ExceptionT::OutOfRange(caller, "(%d,%d) not present", row, col);
	
	/* dummy */
	return ((double*)A->nzval)[0];
}

/* element accessor - READ ONLY */
inline double SuperLU_DISTMatrixT::Element(int row, int col) const
{
	const char caller[] = "SuperLU_DISTMatrixT::Element()";

#if __option(extended_errorcheck)
	/* range checks */
	if (row < 0 || row >= fLocNumEQ) ExceptionT::GeneralFail(caller);
	if (col < 0 || col >= fLocNumEQ) ExceptionT::GeneralFail(caller);
#endif

	NCformat *A = (NCformat*) fA.Store;

	/* look through column col for the given row index */
	int start = A->colptr[col];
	int stop  = A->colptr[col+1];
	for (int i = start; i < stop; ++i)
	{
		if (row == A->rowind[i])
			return ((double*)A->nzval)[i];
	}

	/* otherwise this nonzero wasn't present */
	return 0.0;
}

} /* namespace Tahoe */

#endif /* __TAHOE_MPI__ */
#endif /* __SUPERLU_DIST__ */
#endif /* _SUPER_LU_DIST_MATRIX_T_H_ */
