/* $Id: SolverT.h,v 1.23 2008-05-12 22:32:27 regueiro Exp $ */
/* created: paklein (05/23/1996) */
#ifndef _SOLVER_H_
#define _SOLVER_H_

/* environment */
#include "Environment.h"

/* base class */
#include "iConsoleObjectT.h"
#include "ParameterInterfaceT.h"

/* direct members */
#include "dArrayT.h"
#include "GlobalMatrixT.h"
#include "GlobalT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#ifdef DEM_COUPLING_DEV
#include "FEDEManagerT.h"
#include "FBC_CardT.h"
#endif
#endif

namespace Tahoe {

/* forward declarations */
class FEManagerT;
class iArrayT;
class iArray2DT;
class dMatrixT;
class ElementMatrixT;
template <class TYPE> class RaggedArray2DT;

/** solver base class. This class is responsible for driving the solution
 * procedure over a given series of time increments. The class controls
 * how the solution is determined over a given time step, what determines
 * if the solution has been found, and how to handle failures to find a
 * solution, or other irregularities that occur during the solution procedure.
 * Derived types instantiate different nonlinear solution procedures. */
class SolverT: public iConsoleObjectT, public ParameterInterfaceT
{
public:

	/** global matrix types */
	enum MatrixTypeT {kDiagonalMatrix = 0, /**< diagonal matrix for "matrix-free" methods */
	                   kProfileSolver = 1, /**< symmetric and nonsymmetric profile solvers */
	                      kFullMatrix = 2, /**< full matrix with pivoting */
					           kAztec = 3, /**< sparse, iterative solver */
			                 kSuperLU = 4, /**< sparse, direct solver */
			                 kSPOOLES = 5, /**< sparse, direct solver */
			                 kPSPASES = 6  /**< sparse, direct multi-processor solver */
			                 };

	/** solution status */
	enum SolutionStatusT {kContinue = 0, /**< solution not found after the most recent iteration */
                          kConverged = 1, /**< solution found */
                             kFailed = 2  /**< solution procedure has failed */
                             };

	/** constructors */
	SolverT(FEManagerT& fe_manager, int group);

	/** destructor */
	virtual ~SolverT(void);

	/** factory method. Construct a new instance of a sub-class of SolverT
	 * with the given ParameterInterfaceT name, or return NULL if the name is
	 * not recognized. */
	static SolverT* New(FEManagerT& fe_manager, const char* name, int group);

	/** (re-)configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/* process element group equation data to configure matrix */
	void ReceiveEqns(const iArray2DT& equations) const;
	void ReceiveEqns(const RaggedArray2DT<int>& equations) const;

	/** \name solution steps */
	/*@{*/
	/** start solution step */
	virtual void InitStep(void);

	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int max_iterations) = 0;

#ifdef DEM_COUPLING_DEV
	virtual SolutionStatusT Solve(int max_iterations, FEDEManagerT& fFEDEManager, ArrayT<FBC_CardT>& fGhostFBC) {};
#endif

	/** end solution step */
	virtual void CloseStep(void);

	/** error handler */
	virtual void ResetStep(void);
	
	/** signal time step change. Chance to clear cached values that may depend on the
	 * time increment. */
	virtual void SetTimeStep(double dt);
	/*@}*/

	/** \name assembling the global equation system */
	/*@{*/
	void UnlockRHS(void) { fRHS_lock = kOpen; };
	void LockRHS(void) { fRHS_lock = kLocked; };
	void UnlockLHS(void) { fLHS_lock = kOpen; };
	void LockLHS(void) { fLHS_lock = kLocked; };
	
	void AssembleLHS(const ElementMatrixT& elMat, const nArrayT<int>& eqnos);
	void AssembleLHS(const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
		const nArrayT<int>& col_eqnos);
	void AssembleLHS(const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos);
	void OverWriteLHS(const ElementMatrixT& elMat, const nArrayT<int>& eqnos);
	void DisassembleLHS(dMatrixT& matrix, const nArrayT<int>& eqnos) const;
	void DisassembleLHSDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const;

	void AssembleRHS(const nArrayT<double>& elRes, const nArrayT<int>& eqnos);

	/** assemble forces over the whole system.
	 * \param elRes force vector with length the total number of unknowns */
	void AssembleRHS(const nArrayT<double>& elRes);
	void OverWriteRHS(const dArrayT& elRes, const nArrayT<int>& eqnos);
	void DisassembleRHS(dArrayT& elRes, const nArrayT<int>& eqnos) const;
	/*@}*/

	/* accessor */
	const int& IterationNumber(void) const;

	/** debugging */
	int Check(void) const;
	const dArrayT& RHS(void) const;
	const GlobalMatrixT& LHS(void) const;

	/* return the required equation numbering scope - local by default */
	GlobalT::EquationNumberScopeT EquationNumberScope(void) const;

	/** returns true if solver prefers reordered equations */
	bool RenumberEquations(void);
	
	/** my group */
	int Group(void) const { return fGroup; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** enum for the protected state of the residual and stiffness matrix */
	enum LockStateT {
		  kOpen = 0, /**< open for assembly */
		kLocked = 1, /**< attempts to assemle throw an ExceptionT::kGeneralFail */
		kIgnore = 2  /**< attempts to assemble are silently ignored */
		};

	/** return the magnitude of the residual force */
	double Residual(const dArrayT& force) const;

	/** inner product */	
	double InnerProduct(const dArrayT& v1, const dArrayT& v2) const;

	/** \name method needed for check code GlobalMatrixT::kCheckLHS */
	/*@{*/
	/** return approximate stiffness matrix. Compute and approximate stiffness
	 * matrix by perturbing each degree of freedom in the system. Caller is 
	 * responsible for disposing of the matrix. */
	GlobalMatrixT* ApproximateLHS(const GlobalMatrixT& template_LHS);

	/** compare the two stiffness matricies. Write the results of the comparison
	 * to FEManagerT::Output */
	void CompareLHS(const GlobalMatrixT& ref_LHS, const GlobalMatrixT& test_LHS) const;
	/*@}*/

private:

	/** check matrix type against analysis code, return 1 if
	 * compatible, 0 otherwise */
	int CheckMatrixType(int matrix_type, int analysis_code) const;

	/** set global equation matrix */
	void SetGlobalMatrix(const ParameterListT& params, int check_code);

protected:

	/** the Boss */	
	FEManagerT& fFEManager;

	/** equation group number */
	int fGroup;

	/** \name flags */
	/*@{*/
	int fMatrixType;
	int fPrintEquationNumbers;
	/*@}*/

	/** global equation system */
	/*@{*/
	/** global LHS matrix */
	GlobalMatrixT* fLHS;
	
	/** write protection for the LHS matrix */
	LockStateT fLHS_lock;

	/** runtime flag. Set to true to signal LHS matrix needs to be recalculated. By
	 * default, this is set to true during the call to SolverT::InitStep. */
	bool fLHS_update;
	
	/** perturbation for computing finite difference version of LHS */
	double fPerturbation;

	/** residual */
	dArrayT fRHS;
	
	/** write protection for the RHS vector */
	LockStateT fRHS_lock;
	/*@}*/

	/** runtime data */
	int fNumIteration;
	
	/** eigenvolsver parameters */
	ParameterListT* fEigenSolverParameters;
};

/* inlines */

/* signal time step change */
inline void SolverT::SetTimeStep(double dt)
{
#pragma unused(dt)
}

/* assemble forces over the whole system */
inline void SolverT::AssembleRHS(const nArrayT<double>& elRes)
{
	/* lock state */
	if (fRHS_lock == kIgnore)
		return;
	else if (fRHS_lock == kLocked)
		ExceptionT::GeneralFail("SolverT::AssembleRHS");
	else
		fRHS += elRes;
}

/* assembling the global equation system */
inline void SolverT::AssembleLHS(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
	if (fLHS_lock == kOpen)
		fLHS->Assemble(elMat, eqnos);
	else if (fLHS_lock == kLocked)
		ExceptionT::GeneralFail("SolverT::AssembleLHS", "LHS is locked");
}

inline void SolverT::AssembleLHS(const ElementMatrixT& elMat, const nArrayT<int>& row_eqnos,
	const nArrayT<int>& col_eqnos)
{
	if (fLHS_lock == kOpen)
		fLHS->Assemble(elMat, row_eqnos, col_eqnos);
	else if (fLHS_lock == kLocked)
		ExceptionT::GeneralFail("SolverT::AssembleLHS", "LHS is locked");
}

inline void SolverT::AssembleLHS(const nArrayT<double>& diagonal_elMat, const nArrayT<int>& eqnos)
{
	if (fLHS_lock == kOpen)
		fLHS->Assemble(diagonal_elMat, eqnos);
	else if (fLHS_lock == kLocked)
		ExceptionT::GeneralFail("SolverT::AssembleLHS", "LHS is locked");
}

inline void SolverT::OverWriteLHS(const ElementMatrixT& elMat, const nArrayT<int>& eqnos)
{
	if (fLHS_lock == kOpen)
		fLHS->OverWrite(elMat, eqnos);
	else if (fLHS_lock == kLocked)
		ExceptionT::GeneralFail("SolverT::OverWriteLHS", "LHS is locked");
}

inline void SolverT::DisassembleLHS(dMatrixT& matrix, const nArrayT<int>& eqnos) const
{
	fLHS->Disassemble(matrix, eqnos);
}

inline void SolverT::DisassembleLHSDiagonal(dArrayT& diagonals, const nArrayT<int>& eqnos) const
{
	fLHS->DisassembleDiagonal(diagonals, eqnos);
}

/* debugging */
inline int SolverT::Check(void) const { return fLHS->CheckCode(); }
inline const dArrayT& SolverT::RHS(void) const { return fRHS; }
inline const GlobalMatrixT& SolverT::LHS(void) const {
	if (!fLHS) ExceptionT::GeneralFail("SolverT::LHS", "LHS not set");
	return *fLHS; 
}

/* accessor */
inline const int& SolverT::IterationNumber(void) const { return fNumIteration; }

} /* namespace Tahoe */

#endif /* _SOLVER_H_ */
