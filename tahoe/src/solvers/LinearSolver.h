/* $Id: LinearSolver.h,v 1.4 2002-07-05 22:28:41 paklein Exp $ */
/* created: paklein (05/30/1996) */

#ifndef _LINEAR_SOLVER_H_
#define _LINEAR_SOLVER_H_

/* base class */
#include "SolverT.h"

namespace Tahoe {

/** solver for linear problems */
class LinearSolver: public SolverT
{
public:

	/* constructor */
	LinearSolver(FEManagerT& fe_manager, int group);
	
	/* configure system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);
	
private:

	/* flag to form RHS */
	int fFormLHS;
		// reform conditions:
		// (1) initially
		// (2) if equation system is reconfigured
};

} // namespace Tahoe 
#endif /* _LINEAR_SOLVER_H_ */
