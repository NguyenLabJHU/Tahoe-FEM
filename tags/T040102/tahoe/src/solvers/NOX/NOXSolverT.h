/* $Id: NOXSolverT.h,v 1.1 2002-03-28 16:40:35 paklein Exp $ */
#ifndef _NOX_SOLVER_T_H_
#define _NOX_SOLVER_T_H_

/* base classes */
#include "SolverT.h"
#include "NOXInterfaceT.h"

/* forward declarations */
namespace NOX {
	namespace Solver { 
		class Manager; 
	}
	namespace Tahoe { 
		class Group; 
	}
}

/** interface to the Sandia NOX nonlinear solver library */
class NOXSolverT: public SolverT, public NOXInterfaceT
{
public:

	/** constructor */
	NOXSolverT(FEManagerT& fe_manager);

	/** destructor */
	virtual ~NOXSolverT(void);

	/** generate the solution for the current time sequence */
	virtual void Run(void);

	/** error handler */
	virtual void ResetStep(void);

	/** compute RHS for the given solution vector x.
	 * \param x solution vector to apply the system
	 * \param rhs returns with the residual associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeRHS(const dArrayT& x, dArrayT& rhs);
  
	/** compute the Jacobian given the specified solution vector x.  
	 * \param x solution vector to apply the system
	 * \param jacobian returns with the Jacobian associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeJacobian(const dArrayT& x, GlobalMatrixT& jacobian);

private:

	/* iteration status flags */
	enum SolutionStatusT {kContinue = 0,
                         kConverged = 1,
                            kFailed = 2};

	/** handle situation if solution for the current time increment
	 * is successfully determined.
	 * \return true if step is ready to be closed */
	SolutionStatusT DoConverged(void);

	/** handle situation if solution for the current time increment
	 * could not be determined */
	void DoNotConverged(void);

	/** divert output for iterations */
	void InitIterationOutput(void);

	/** restore normal output */
	void CloseIterationOutput(void);

protected:
	
	/** NOX solver */
	NOX::Solver::Manager* fNOX;

	/** seed NOX group sent to the solver */
	NOX::Tahoe::Group* fGroup;
	
private:

	/* parameters */
	int fIterationOutputIncrement; /**< "movies" of convergence steps */

	/* runtime data */
	GlobalMatrixT* fSeed_LHS;  /**< seed LHS matrix */
	int fIterationOutputCount; /**< output count for NOXSolverT::fIterationOutputIncrement */
};

#endif /* _NOX_SOLVER_T_H_ */
