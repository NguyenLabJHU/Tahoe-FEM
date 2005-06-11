/* $Id: PenaltySphereT.h,v 1.8 2005-02-22 00:10:19 rjones Exp $ */
/* created: paklein (04/30/1998) */
#ifndef _PENATLY_SPHERE_T_H_
#define _PENATLY_SPHERE_T_H_

/* base class */
#include "PenaltyRegionT.h"

/* direct members */
#include "ElementMatrixT.h"

namespace Tahoe {


/** spherical rigid barrier enforced with a penalized constraint */
class PenaltySphereT: public PenaltyRegionT
{
public:

	/* constructor */
	PenaltySphereT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

protected:
	/** nodal areas */
	double fRadius;

	/** \name workspace */
	/*@{*/
	dArrayT        fv_OP; /**< vector from center to contact node */
	ElementMatrixT fLHS;  /**< tangent matrix */
	dArrayT        fd_sh; //shallow
	iArrayT        fi_sh; //shallow
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PENATLY_SPHERE_T_H_ */
