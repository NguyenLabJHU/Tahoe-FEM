/* $Id: NOXSolverT.h,v 1.5 2002-07-03 23:11:11 paklein Exp $ */
#ifndef _NOX_SOLVER_T_H_
#define _NOX_SOLVER_T_H_

/* optional */
#ifdef __NOX__

/* base classes */
#include "SolverT.h"
#include "SolverInterfaceT.h"

/* forward declarations */
namespace NOX {
	namespace Parameter {
		class List;
	}
	namespace Solver {
		class Manager;
	}
	namespace Status {
		class AbsResid;
//		class RelResid;
		class MaxIters;
	}
} // namespace NOX

namespace Tahoe {

/** interface to the Sandia NOX nonlinear solver library. */
class NOXSolverT: public SolverT, protected SolverInterfaceT
{
public:

	/** constructor 
	 * \param unknowns_order the time derivative of the field that is
	 *        treated as the unknown field by the solver */
	NOXSolverT(FEManagerT& fe_manager, int group, int unknowns_order);

	/** destructor */
	virtual ~NOXSolverT(void);

	/** solve the system over the current time increment.
	 * \param num_iterations maximum number of iterations to execute. Hitting this limit
	 *        does not signal a SolverT::kFailed status, unless solver's internal parameters
	 *        also indicate the solution procedure has failed.
	 * \return one of SolverT::IterationsStatusT */
	virtual SolutionStatusT Solve(int num_iterations);

	/** error handler */
	virtual void ResetStep(void);

	/** (re-)configure the global equation system */
	virtual void Initialize(int tot_num_eq, int loc_num_eq, int start_eq);

protected:
	
	/*@{ \name concrete implementation of the SolverInterfaceT.*/
	/** compute RHS for the given update to the solution. 
	 * \param x the solution vector to apply the system
	 * \param rhs returns with the residual associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeRHS(const dArrayT& x, dArrayT& rhs);
  
	/** compute the Jacobian associated with the solution used for the most
	 * recent call to NOXInterfaceT::computeRHS.  
	 * \param jacobian returns with the Jacobian associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeJacobian(GlobalMatrixT& jacobian);
	/*@}*/

private:

	/** run solver */                            
	SolutionStatusT Solve(NOX::Solver::Manager& nox);

	/** divert output for iterations */
	void InitIterationOutput(void);

	/** restore normal output */
	void CloseIterationOutput(void);

protected:

	/** parameter list for the NOX solver. The parameters may specify the method:
	 *  <ul>
	 *  <li> "Nonlinear %Solver" - Name of the solver method. Valid choices are
	 *  <ul> 
	 *  <li> "%Newton" (NOX::Solver::Newton) [Default]
	 *  <li> "Nonlinear CG" (NOX::Solver::NonlinearCG)
	 *  <li> "Trust Region" (NOX::Solver::TrustRegion)
	 *  </ul> </ul> 
	 *
	 * Additional parameters can be added to the sublists "Direction" and
	 * "Line Search", see NOX::Direction::Manager and NOX::Linesearch::Manager
	 * for more information. */
	NOX::Parameter::List* fNOXParameters;
	
private:

	/** order of the field solved by the solver */
	int fUnknownsOrder;

	/** \name stopping criteria */
	/*@{*/
	int    fMaxIterations; /**< maximum number of iterations */
	double fAbsResidual;   /**< magnitude of the residual */
	double fRelResidual;   /**< relative magnitude of the residual */
	/*@}*/

	/** \name step growth parameters. 	
	 * Step cuts occur if the solution fails to converge. Step growth is controlled
	 * by a dependent, two-parameter condition. */
	/*@{*/
	int fQuickSolveTol;  /**< iterations considered "easy" solution */
	int fQuickSeriesTol; /**< number "easy" solutions before step increase */
	/*@}*/

	/** "movies" of convergence steps */
	int fIterationOutputIncrement; 

	/** \name runtime data */
	/*@{*/
	int fIterationOutputCount; /**< output count for NOXSolverT::fIterationOutputIncrement */
	dArrayT fLastSolution;     /**< the last solution applied to the system */
	/*@}*/
};

} // namespace Tahoe
#endif /* __NOX__ */
#endif /* _NOX_SOLVER_T_H_ */
