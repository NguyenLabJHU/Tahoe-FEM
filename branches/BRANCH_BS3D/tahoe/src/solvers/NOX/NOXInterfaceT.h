/* $Id: NOXInterfaceT.h,v 1.2 2002-07-02 19:57:16 cjkimme Exp $ */
#ifndef NOX_TAHOE_INTERFACE_H
#define NOX_TAHOE_INTERFACE_H
  
/* forward declarations */

namespace Tahoe {

class dArrayT;
class GlobalMatrixT;  
  
/** abstract interface for allow NOX to get residual and Jacobian
 * information from Tahoe.
 * \note Due to the conventions in Tahoe, NOXInterfaceT::computeRHS must be called 
 *       before the corresponding call to NOXInterfaceT::computeJacobian. This is
 *       enforced by not allowing for the new solution vector when asking for the
 *       Jacobian. */
class NOXInterfaceT {

  public:
 
  	/** constructor */
	NOXInterfaceT(void) { };

	/** destructor */
	virtual ~NOXInterfaceT(void) { };

	/** compute RHS for the given solution vector x. 
	 * \param x solution vector to apply the system
	 * \param rhs returns with the residual associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeRHS(const dArrayT& x, dArrayT& rhs) = 0;
  
	/** compute the Jacobian associated with the solution used for the most
	 * recent call to NOXInterfaceT::computeRHS.  
	 * \param jacobian returns with the Jacobian associated with the given solution vector
	 * \return true if computation was successful */
	virtual bool computeJacobian(GlobalMatrixT& jacobian) = 0;
};

} // namespace Tahoe 
#endif /* NOX_TAHOE_INTERFACE_H */
