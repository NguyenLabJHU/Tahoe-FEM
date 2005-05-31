/* $Id: LinearSolver.h,v 1.8.20.2 2005-05-31 06:08:55 paklein Exp $ */
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

	/** constructor */
	LinearSolver(FEManagerT& fe_manager, int group);
	
	/** configure system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

	/** start solution step */
	virtual void InitStep(void);
	
	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);

	/** signal time step change. Chance to clear cached values that may depend on the
	 * time increment. LinearSolver::SetTimeStep triggers recalculation of the LHS
	 * matrix because some time integrators use an effective mass matrix that is a function
	 * of the time increment. */
	virtual void SetTimeStep(double dt);
	
private:

	/* flag to form RHS */
	int fFormLHS;
		// reform conditions:
		// (1) initially
		// (2) if equation system is reconfigured
		// (3) when the time step changes
};

} // namespace Tahoe 
#endif /* _LINEAR_SOLVER_H_ */
