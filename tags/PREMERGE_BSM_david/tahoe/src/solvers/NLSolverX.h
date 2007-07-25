/* $Id: NLSolverX.h,v 1.8 2004-09-09 23:54:55 paklein Exp $ */
/* created: paklein (08/25/1996) */
#ifndef _NL_SOLVER_X_H_
#define _NL_SOLVER_X_H_

/* base class */
#include "NLSolver.h"

namespace Tahoe {

/* forward declarations */
class CCSMatrixT;
class CCNSMatrixT;

/** nonlinear solver methods testbed */
class NLSolverX: public NLSolver
{
public:

	/** constructor */
	NLSolverX(FEManagerT& fe_manager, int group);

	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);

	/* error handler */
	virtual void ResetStep(void);

protected:

#if 0
	/* relax system - reform tangent at newtancount intervals */
	virtual SolutionStatusT Relax(int newtancount = 1); 
#endif

	/** things to do if the solver converges */
	SolutionStatusT DoConverged(void);

private:

	/* parameters */
	int fMaxNewTangents;  // max number of new tangents per step
	int fMaxTangentReuse; // max number iterations with same tangent
	int fMinFreshTangents;// number of new tan's at start of step (-1 for auto)
	int fCheckNegPivots; // automatic new tangent for negative pivots
	double fTrustTol;     // conv below which an extra 3 new tangents is allowed
	double fMinResRatio;  // min quasi-Newton reduction ratio: r_i+1/r_i

	/* runtime parameters */
	int fFormNewTangent; // 1 => new tangent, 0 => re-use
	int fNumNewTangent;

	/* store last update for recovery */
	dArrayT fLastUpdate;
	
	/* dynamically casted */
	CCSMatrixT*  pCCS;
	CCNSMatrixT* pCCNS;
};

} // namespace Tahoe 
#endif /* _NL_SOLVER_X_H_ */
