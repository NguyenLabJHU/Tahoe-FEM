/* $Id: NLSolver.h,v 1.2.2.3 2002-04-30 08:22:05 paklein Exp $ */
/* created: paklein (07/09/1996) */

#ifndef _NL_SOLVER_H_
#define _NL_SOLVER_H_

/* base class */
#include "SolverT.h"

/** nonlinear Newton solver. */
class NLSolver: public SolverT
{
public:

	/* constructor */
	NLSolver(FEManagerT& fe_manager, int group);
	
	/** solve the system over the current time increment */
	virtual int Solve(void);	

	/* error handler */
	virtual void ResetStep(void);

#ifdef _MSC_VER
	/* iteration status flags */
	enum IterationStatusT {kContinue = 0,
                          kConverged = 1,
                             kFailed = 2};
#endif

protected:

#ifndef _MSC_VER
	/* iteration status flags */
	enum IterationStatusT {kContinue = 0,
                          kConverged = 1,
                             kFailed = 2};
#endif

	/* apply system update (socket for line searching), pass NULL
	 * for residual if not available */
	virtual void Update(const dArrayT& update, const dArrayT* residual);

	/* relax system - reform tangent at newtancount intervals */
	virtual IterationStatusT Relax(int newtancount = 1);

	/* advance to next load step. Returns 0 if there are no more
	 * steps. Overload to add class dependent initializations. */
//	virtual int Step(void);

	/* returns the appropriate iteration status flag for
	 * the given error measurement, based on the current
	 * iteration number, convergence tolerance, etc. */
	IterationStatusT ExitIteration(double error);

	/* form and solve the equation system - returns the magnitude of the
	 * residual */
	virtual double SolveAndForm(bool newtangent);

	/* divert output for iterations */
	void InitIterationOutput(void);
	void CloseIterationOutput(void);

protected:

	/** things to do if the solver converges */
	IterationStatusT DoConverged(void);

protected:

	/* error management parameters */	
	int    fMaxIterations;  // maximum number of iterations per step
	double fZeroTolerance;  // absolute convergence tolerance
	double fTolerance;		// relative convergence tolerance
	double fDivTolerance;   // tolerance for a diverging solution
	int    fQuickSolveTol;  // iterations considered "easy" solution
	int    fQuickSeriesTol; // "easy" solutions before step increase
	int    fIterationOutputIncrement; // "movies" of convergence steps

	/* runtime error management data */
	double fError0;
	int	   fQuickConvCount;
	int    fIterationOutputCount;

	/* output control */
	int fVerbose;
};

#endif /* _NL_SOLVER_H_ */
