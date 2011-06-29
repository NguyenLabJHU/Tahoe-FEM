/* $Id: SuperLU_DISTMatrixT.h,v 1.8 2005-08-01 03:26:30 paklein Exp $ */
#ifndef _SUPER_LU_DIST_MATRIX_T_H_
#define _SUPER_LU_DIST_MATRIX_T_H_

/* library support */
#ifdef __SUPERLU_DIST__
#ifdef __TAHOE_MPI__

/* base class */
#include "GlobalMatrixT.h"

/* SuperLU type definitions */
#include "superlu_ddefs.h"

/* direct members */
#include "AutoArrayT.h"

namespace Tahoe {

/* forward declarations */
class MSRBuilderT;

/** interface to SuperLU 2.0 parallel linear solver */
class SuperLU_DISTMatrixT: public GlobalMatrixT
{
public:

	/** constructor */
	SuperLU_DISTMatrixT(ostream& out, int check_code, bool print_stat, 
		IterRefine_t refine, const CommunicatorT& comm);

	/** copy constructor */
	SuperLU_DISTMatrixT(const SuperLU_DISTMatrixT& rhs);

	/** destructor */
	~SuperLU_DISTMatrixT(void);

	/** enum conversion */
	static IterRefine_t int2IterRefine_t(int i);

	/** set the internal matrix structure */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** write information to output stream after SuperLU_DISTMatrixT::Initialize
	 * has been called */
	virtual void Info(ostream& out);
	
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

	/* number scope and reordering */
	virtual EquationNumberScopeT EquationNumberScope(void) const;
	virtual bool RenumberEquations(void) const;

	/** return the form of the matrix */
	virtual GlobalT::SystemTypeT MatrixType(void) const { return GlobalT::kNonSymmetric; };

	/** assignment operator */
	SuperLU_DISTMatrixT& operator=(const SuperLU_DISTMatrixT& rhs);

	/** return a clone of self */
	virtual GlobalMatrixT* Clone(void) const;

protected:

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

	/** element accessor. Returns a pointer to the given element in the matrix
	 * or NULL if the element is not in the set of non-zero elements. The search
	 * for columns is performed through recursive bisection of the indicies of
	 * columns with non-zero elements in the specified row. */
	double* operator()(int row, int col);

private:

	/** \name clean up methods */
	/*@{*/
	void FreeLUstruct(int dim, gridinfo_t& grid, LUstruct_t& lu_struct) const;
	void FreeScalePermstruct(ScalePermstruct_t& scale_struct) const;
	/*@}*/

protected:

	/** matrix structure builder */
	MSRBuilderT* fBuilder;

	/** \name matrix and factors in SuperLU formats */
	/*@{*/
	SuperMatrix fA;
	AutoArrayT<int> frowptr;
	AutoArrayT<int> fcolind;
	AutoArrayT<double> fnzval;

    ScalePermstruct_t fScalePermstruct;
    LUstruct_t fLUstruct;
    SOLVEstruct_t fSOLVEstruct;

	/** copy of SuperLU_DISTMatrixT::fcolind */
	AutoArrayT<int> fcolind2;

	/** copy of SuperLU_DISTMatrixT::frowptr */
	AutoArrayT<int> frowptr2;
	/*@}*/

	/** SuperLU_DIST options */
	superlu_options_t foptions;
	
	/** process grid information */
	gridinfo_t fgrid;

	/** \name factorization flags */
	/*@{*/
	/** true if matrix has been symbolically factorized */
	bool fIsSymFactorized;

	/** true if matrix has been numerically factorized */
	bool fIsNumFactorized;
	/*@}*/
};

/* element accessor */
inline double* SuperLU_DISTMatrixT::operator()(int row, int col)
{
	const char caller[] = "SuperLU_DISTMatrixT::operator()";
	int loc_row = row - fStartEQ + 1; /* fStartEQ is 1,... */

#if __option(extended_errorcheck)
	/* range checks */
	if (loc_row < 0 || loc_row >= fLocNumEQ) ExceptionT::OutOfRange(caller, "row_loc %d < 0 or >= %d", loc_row, fLocNumEQ);
	if (col < 0 || col >= fTotNumEQ) ExceptionT::OutOfRange(caller, "col %d < 0 || >= %d", col, fTotNumEQ);
#endif

	/* non-zero columns */
	int r_dex = frowptr[loc_row];
	int* pcolind = fcolind.Pointer(r_dex);

	/* range */
	int min = 0;
	int max = frowptr[loc_row+1] - r_dex - 1;

	/* is last value */
	if (pcolind[max] == col)
		return fnzval.Pointer(r_dex + max);

	/* bisection */
	int c_dex = (max + min)/2;
	int c_hit = pcolind[c_dex];
	while (c_hit != col && max != min+1) {
	
		/* shift bounds */
		if (c_hit > col)
			max = c_dex;
		else /* painds[cdex] < col */
			min = c_dex;
	
		/* bisect */
		c_dex = (max + min)/2;
		c_hit = pcolind[c_dex];
	}
	
	/* found */
	if (c_hit == col)
		return fnzval.Pointer(r_dex + c_dex);
	else /* not found */
		return NULL;
}

/* enum conversion */
inline IterRefine_t SuperLU_DISTMatrixT::int2IterRefine_t(int i)
{
	switch (i)
	{
		case NOREFINE:
			return NOREFINE;
		case SINGLE:
			return SINGLE;
		case DOUBLE:
			return DOUBLE;
		case EXTRA:
			return EXTRA;
		default:
			ExceptionT::GeneralFail("SuperLU_DISTMatrixT::int2IterRefine_t", 
				"unrecognized %d", i);
	}
	return NOREFINE;
}

} /* namespace Tahoe */

#endif /* __TAHOE_MPI__ */
#endif /* __SUPERLU_DIST__ */
#endif /* _SUPER_LU_DIST_MATRIX_T_H_ */
