/* $Id: PenaltyCylinderT.h,v 1.3 2004-07-15 08:31:15 paklein Exp $ */
#ifndef _PENALTY_CYLINDER_T_H_
#define _PENALTY_CYLINDER_T_H_

/* base class */
#include "PenaltyRegionT.h"

/* direct members */
#include "ElementMatrixT.h"

namespace Tahoe {

/** rigid cylinder enforced with a penalty formulation */
class PenaltyCylinderT: public PenaltyRegionT
{
public:

	/** constructor */
	PenaltyCylinderT(void);

	/** form of tangent matrix */
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

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/
	
protected:

	/** accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce);

protected:

	/** sphere radius */
	double fRadius;
	
	/** cylinder direction */
	dArrayT fDirection;
	
	/** \name workspace */
	/*@{*/
	dArrayT        fv_OP; /**< vector from center to contact node */
	dArrayT        fR;    /**< vector normal to the cylinder */
	ElementMatrixT fLHS;  /**< tangent matrix */
	dArrayT        fd_sh; /**< shallow */
	iArrayT        fi_sh; /**< shallow */
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PENALTY_CYLINDER_T_H_ */
