/* $Id: NLSolver.h,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (07/09/1996)                                          */

#ifndef _NL_SOLVER_H_
#define _NL_SOLVER_H_

/* base class */
#include "SolverT.h"

class NLSolver: public SolverT
{
public:

	/* constructor */
	NLSolver(FEManagerT& fe_manager);
	
	/* generate the solution for the current time sequence */
	 virtual void Run(void);

	/* error handler */
	virtual void ResetStep(void);

protected:

	/* iteration status flags */
	enum IterationStatusT {kContinue = 0,
                          kConverged = 1,
                             kFailed = 2};

	/* apply system update (socket for line searching), pass NULL
	 * for residual if not available */
	virtual void Update(const dArrayT& update, const dArrayT* residual);

	/* relax system - reform tangent at newtancount intervals */
	virtual IterationStatusT Relax(int newtancount = 1);

	/* advance to next load step. Returns 0 if there are no more
	 * steps. Overload to add class dependent initializations. */
	virtual int Step(void);

	/* returns the appropriate iteration status flag for
	 * the given error measurement, based on the current
	 * iteration number, convergence tolerance, etc. */
	IterationStatusT ExitIteration(double error);

	/* form and solve the equation system - returns the magnitude of the
	 * residual */
	virtual double SolveAndForm(bool newtangent);
	
	/* handlers */
	virtual IterationStatusT DoConverged(void);
	virtual void DoNotConverged(void);

	/* divert output for iterations */
	void InitIterationOutput(void);
	void CloseIterationOutput(void);

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
