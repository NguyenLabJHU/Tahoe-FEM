/* $Id: SolverT.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (05/23/1996)                                          */

#ifndef _SOLVER_H_
#define _SOLVER_H_

/* environment */
#include "Environment.h"

/* base class */
#include "iConsoleObjectT.h"

/* direct members */
#include "dArrayT.h"
#include "GlobalMatrixT.h" // for the inlines
#include "GlobalT.h"

/* forward declarations */
class FEManagerT;
class iArrayT;
class iArray2DT;
class dMatrixT;
class ElementMatrixT;
template <class TYPE> class RaggedArray2DT;

class SolverT: public iConsoleObjectT
{
public:

	/* nonlinear solver codes */
	enum SolverTypeT {kNewtonSolver = 0, // standard Newton solver
                   kK0_NewtonSolver = 1, // initial tangent, Newton solver
                   kModNewtonSolver = 2, // modified Newton solver (development)
                    kExpCD_DRSolver = 3, // central difference, dynamic relaxation
                   kNewtonSolver_LS = 4, // Newton solver with line search
                      kPCGSolver_LS = 5, // preconditioned, nonlinear conjugate gradient
                  kiNewtonSolver_LS = 6};// interactive Newton solver (with line search)

	/* global matrix types */
	enum MatrixTypeT {kDiagonalMatrix = 0,
	                   kProfileSolver = 1,  // symmetric and nonsymmetric
	                      kFullMatrix = 2,  // with pivoting
					           kAztec = 3,  // sparse, iterative solver
			            kSparseDirect = 4,  // sparse, direct solver: SuperLU
			                 kSPOOLES = 5}; // sparse, direct solver: symbolic factorization

	/* constructor */
	SolverT(FEManagerT& fe_manager);

	/* destructor */
	virtual ~SolverT(void);

	/* (re-)configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/* process element group equation data to configure matrix */
	void ReceiveEqns(const iArray2DT& equations) const;
	void ReceiveEqns(const RaggedArray2DT<int>& equations) const;

	/* Generate the solution for the current time sequence */
	virtual void Run(void) = 0;

	/* error handler */
	virtual void ResetStep(void);
	
	/* assembling the global equation system */
	void AssembleLHS(const ElementMatrixT& elMat, const iArrayT& eqnos);
	void OverWriteLHS(const ElementMatrixT& elMat, const iArrayT& eqnos);
	void DisassembleLHS(dMatrixT& elMat, const iArrayT& eqnos) const;

	void AssembleRHS(const dArrayT& elRes, const iArrayT& eqnos);
	void OverWriteRHS(const dArrayT& elRes, const iArrayT& eqnos);
	void DisassembleRHS(dArrayT& elRes, const iArrayT& eqnos) const;

	/* accessor */
	const int& IterationNumber(void) const;

	/* debugging */
	int Check(void) const;
	const dArrayT& RHS(void) const;

	/* return the required equation numbering scope - local by default */
	GlobalT::EquationNumberScopeT EquationNumberScope(void) const;

	/* returns true if solver prefers reordered equations */
	bool RenumberEquations(void);

protected:

	/* advance to next load step. Returns 0 if there are no more
	 * steps. Overload to add class dependent initializations */
	virtual int Step(void);

	/* return the magnitude of the residual force */
	double Residual(const dArrayT& force) const;

	/* inner product */	
	double InnerProduct(const dArrayT& v1, const dArrayT& v2) const;

private:

	/* check matrix type against analysis code, return 1 if
	 * compatible, 0 otherwise */
	int CheckMatrixType(int matrix_type, int analysis_code) const;

	/* set global equation matrix */
	void SetGlobalMatrix(int matrix_type, int check_code);
		 	
protected:
	
	FEManagerT& fFEManager;

	/* flags */
	int fMatrixType;
	int fPrintEquationNumbers;
	
	/* global equation system */
	GlobalMatrixT* fLHS;	
	dArrayT        fRHS;

	/* runtime data */
	int fNumIteration;
};

/* inlines */

/* assembling the global equation system */
inline void SolverT::AssembleLHS(const ElementMatrixT& elMat, const iArrayT& eqnos)
{
	fLHS->Assemble(elMat, eqnos);
}

inline void SolverT::OverWriteLHS(const ElementMatrixT& elMat, const iArrayT& eqnos)
{
	fLHS->OverWrite(elMat, eqnos);
}

inline void SolverT::DisassembleLHS(dMatrixT& elMat, const iArrayT& eqnos) const
{
	fLHS->Disassemble(elMat, eqnos);
}

/* debugging */
inline int SolverT::Check(void) const { return fLHS->CheckCode(); }
inline const dArrayT& SolverT::RHS(void) const { return fRHS; }

/* accessor */
inline const int& SolverT::IterationNumber(void) const { return fNumIteration; }

#endif /* _SOLVER_H_ */
