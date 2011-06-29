/* $Id: SolverInterfaceT.h,v 1.3 2002-07-05 22:28:41 paklein Exp $ */
#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H
  
namespace Tahoe {

/* forward declarations */
class dArrayT;
class GlobalMatrixT;  
  
/** abstract interface to allow nonlinear solvers to get residual and Jacobian
 * information from Tahoe.
 * \note Due to the conventions in Tahoe, SolverInterfaceT::computeRHS must be called 
 *       before the corresponding call to SolverInterfaceT::computeJacobian. This is
 *       enforced by not allowing for the new solution vector when asking for the
 *       Jacobian. */
class SolverInterfaceT {

  public:
 
  	/** constructor */
	SolverInterfaceT(void) { };

	/** destructor */
	virtual ~SolverInterfaceT(void) { };

	/** compute RHS for the given update to the solution. 
	 * \param x the solution vector to apply the system
	 * \param rhs returns with the residual associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeRHS(const dArrayT& x, dArrayT& rhs) = 0;
  
	/** compute the Jacobian associated with the solution used for the most
	 * recent call to SolverInterfaceT::computeRHS.  
	 * \param jacobian returns with the Jacobian associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeJacobian(GlobalMatrixT& jacobian) = 0;

};

} // namespace Tahoe 
#endif /* SOLVER_INTERFACE_H */
