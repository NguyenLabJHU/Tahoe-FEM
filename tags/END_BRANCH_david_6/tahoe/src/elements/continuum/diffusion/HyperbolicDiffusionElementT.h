/* $Id: HyperbolicDiffusionElementT.h,v 1.1 2004-10-20 21:44:42 paklein Exp $ */
#ifndef _HYPERBOLIC_DIFFUSE_T_H_
#define _HYPERBOLIC_DIFFUSE_T_H_

/* base class */
#include "DiffusionElementT.h"

namespace Tahoe {

/** hyperbolic, linear diffusion element governed by the balance equation
 \f[
 	c \tau T,_{tt} + c T,_{t} = \mathbf{k} \nabla^2 T + r
 \f]
 * where heat flux and temperature are related by the Maxwell-Cattaneo relation
 \f[
	\tau \mathbf{q},_{t} + \mathbf{q} = - \mathbf{k} \nabla T
 \f]
  **/
class HyperbolicDiffusionElementT: public DiffusionElementT
{
public:

	/** constructor */
	HyperbolicDiffusionElementT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** allocate and initialize local arrays */
	virtual void SetLocalArrays(void);

	/** construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/** form the residual force vector */
	virtual void RHSDriver(void);

protected:

	/** relaxation time */
	double fTau;

	/** temperature "acceleration" */
	LocalArrayT fLocAcc;
};

} /* namespace Tahoe */

#endif /* _HYPERBOLIC_DIFFUSE_T_H_ */
