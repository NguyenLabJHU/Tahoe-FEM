/* $Id: DRSolver.h,v 1.5 2004-07-15 08:31:50 paklein Exp $ */
/* created: PAK/CBH (10/03/1996) */

#ifndef _DRSOLVER_H_
#define _DRSOLVER_H_

/* base class */
#include "NLSolver.h"

namespace Tahoe {

/* forward declarations */
class CCSMatrixT;

/** nonlinear solver using dynamic relaxation */
class DRSolver: public NLSolver
{
public:

	/** constructor */
	DRSolver(FEManagerT& fe_manager, int group);

	/* configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);
	
	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);
	
private:

	/* compute the pseudo-mass */
	void ComputeMass(void);
	void ComputeVelocity(void);
	void ComputeDamping(void);

private:

	dArrayT fMass;
	dArrayT fVel;
	dArrayT fDisp;
	double  fDamp;
	
	/* work space */
	dArrayT fKd;		

	int fNumEquations;
	CCSMatrixT* fCCSLHS; //dynamically casted base class data
	 	  	
};

} // namespace Tahoe 
#endif /* _DRSOLVER_H_ */
