/* $Id: NLSolver.h,v 1.10 2004-07-15 08:31:50 paklein Exp $ */
/* created: paklein (07/09/1996) */

#ifndef _NL_SOLVER_H_
#define _NL_SOLVER_H_

/* base class */
#include "SolverT.h"

namespace Tahoe {

/** nonlinear Newton solver. */
class NLSolver: public SolverT
{
public:

	/** constructors */
	NLSolver(FEManagerT& fe_manager, int group);
	
	/** \name solution steps */
	/*@{*/
	/** start solution step */
	virtual void InitStep(void);

	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int max_iterations);

	/** end solution step */
	virtual void CloseStep(void);

	/** error handler */
	virtual void ResetStep(void);	
	/*@}*/

	/** (re-)set the reference error */
	void SetReferenceError(double error);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* apply system update (socket for line searching), pass NULL
	 * for residual if not available */
	virtual void Update(const dArrayT& update, const dArrayT* residual);

	/* relax system - reform tangent at newtancount intervals */
	virtual SolutionStatusT Relax(int newtancount = 1);

	/** returns the appropriate iteration status flag for
	 * the given error measurement, based on the current
	 * iteration number, convergence tolerance, etc. */
	SolutionStatusT ExitIteration(double error, int iteration);

	/* form and solve the equation system - returns the magnitude of the
	 * residual */
	virtual double SolveAndForm(int& iteration);

	/* divert output for iterations */
	void InitIterationOutput(void);
	void CloseIterationOutput(void);

protected:

	/** things to do if the solver converges */
	SolutionStatusT DoConverged(void);

protected:

	/** \name error management parameters */	
	/*@{*/
	int    fMaxIterations;  /**< maximum number of iterations per step */
	int    fMinIterations;  /**< minimum number of iterations per step */
	double fZeroTolerance;  /**< absolute convergence tolerance */
	double fTolerance;		/**< relative convergence tolerance */
	double fDivTolerance;   /**< tolerance for a diverging solution */
	int    fQuickSolveTol;  /**< iterations considered "easy" solution */
	int    fQuickSeriesTol; /**< "easy" solutions before step increase */
	int    fIterationOutputIncrement; /**< "movies" of convergence steps */
	/*@}*/

	/** \name runtime error management data */
	/*@{*/
	double fError0;
	int	   fQuickConvCount;
	int    fIterationOutputCount;
	/*@}*/

	/* output control */
	int fVerbose;
};

} // namespace Tahoe 
#endif /* _NL_SOLVER_H_ */
