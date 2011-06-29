/* $Id: PenaltyContact2DT.h,v 1.7 2004-07-15 08:26:08 paklein Exp $ */
/* created: paklein (12/11/1997) */
#ifndef _PENALTY_CONTACT2D_T_H_
#define _PENALTY_CONTACT2D_T_H_

/* base classes */
#include "Contact2DT.h"

namespace Tahoe {

/** 2D penalty contact element */
class PenaltyContact2DT: public Contact2DT
{
public:

	/** constructor */
	PenaltyContact2DT(const ElementSupportT& support);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** construct the residual force vector */
	virtual void RHSDriver(void);
		
protected:

	/** penalty "stiffness" */
	double fK;

	/** \name element coords and displacements */
	/*@{*/
	dArray2DT fElCoord;
	dArray2DT fElDisp;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PENALTY_CONTACT2D_T_H_ */
