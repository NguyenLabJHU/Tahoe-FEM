/* $Id: NLK0Solver.h,v 1.3 2002-06-08 20:20:55 paklein Exp $ */
/* created: paklein (10/01/1996) */

#ifndef _NL_K0_SOLVER_H_
#define _NL_K0_SOLVER_H_

/* base class */
#include "NLSolver.h"

/* direct members */
#include "CCSMatrixT.h"

/** solver that stores and re-uses the last stiffness matrix
 * which was positive definite. The matrix type must be CCSMatrixT. */
class NLK0Solver: public NLSolver
{
public:

	/** constructor */
	NLK0Solver(FEManagerT& fe_manager, int group);

protected:

	/** form and solve - returns the magnitude of the residual */
	virtual double SolveAndForm(bool junk);

private:

	/* casted pointer to global matrix */
	CCSMatrixT*	pCCSLHS;
	
	/* last positive definite tangent (factorized) */
	int			fFormTangent; //no switching back, yet
	CCSMatrixT	fLastTangent;

	/* full update vector */
	dArrayT		fUpdate;

};

#endif /* _NL_K0_SOLVER_H_ */
