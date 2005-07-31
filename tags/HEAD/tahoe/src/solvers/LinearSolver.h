/* $Id: LinearSolver.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (05/30/1996)                                          */

#ifndef _LINEAR_SOLVER_H_
#define _LINEAR_SOLVER_H_

/* base class */
#include "SolverT.h"

class LinearSolver: public SolverT
{
public:

	/* constructor */
	LinearSolver(FEManagerT& fe_manager);
	
	/* configure system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/* solve for the current time sequence */
	 virtual void Run(void);
	
private:

	/* flag to form RHS */
	int fFormLHS;
		// reform conditions:
		// (1) initially
		// (2) if equation system is reconfigured
};

#endif /* _LINEAR_SOLVER_H_ */
