/* $Id: iNLSolver_LS.h,v 1.5 2002-07-05 22:28:41 paklein Exp $ */
/* created: paklein (01/01/2001) */

#ifndef _I_NL_SOLVER_LS_H_
#define _I_NL_SOLVER_LS_H_

/* base class */
#include "NLSolver_LS.h"

/* direct members */
#include "dArray2DT.h"

namespace Tahoe {

/** nonlinear Newton solver with interactive console */
class iNLSolver_LS: public NLSolver_LS
{
public:

	/* constructor */
	iNLSolver_LS(FEManagerT& fe_manager, int group);

	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);
	
	/* execute commands */
	virtual bool iDoCommand(const CommandSpecT& command, StringT& line);

private:

	/* apply flags */
	virtual void Update(const dArrayT& update, const dArrayT* residual);

	/* commands */
	bool DoStep(int max_steps);
	bool DoInitStep(void);
	SolutionStatusT DoIterate(int max_count);
	
private:

	/* ON/OFF flags */
	bool fFormTangent;	
	bool fLineSearch;	

	/* solution states */
	SolutionStatusT fIterationStatus;
};

} // namespace Tahoe 
#endif /* _I_NL_SOLVER_LS_H_ */
