/* $Id: NLK0Solver.h,v 1.6 2002-12-13 02:42:55 paklein Exp $ */
/* created: paklein (10/01/1996) */

#ifndef _NL_K0_SOLVER_H_
#define _NL_K0_SOLVER_H_

/* base class */
#include "NLSolver.h"

/* direct members */
#include "CCSMatrixT.h"

namespace Tahoe {

/** solver that stores and re-uses the last stiffness matrix
 * which was positive definite. The matrix type must be CCSMatrixT. */
class NLK0Solver: public NLSolver
{
public:

	/** constructor */
	NLK0Solver(FEManagerT& fe_manager, int group);

protected:

	/** form and solve - returns the magnitude of the residual */
	virtual double SolveAndForm(bool junk, bool clear_LHS);

private:

	/* casted pointer to global matrix */
	CCSMatrixT*	pCCSLHS;
	
	/* last positive definite tangent (factorized) */
	int			fFormTangent; //no switching back, yet
	CCSMatrixT	fLastTangent;

	/* full update vector */
	dArrayT		fUpdate;

};

} // namespace Tahoe 
#endif /* _NL_K0_SOLVER_H_ */
