/* $Id: PenaltyContactDrag2DT.h,v 1.2 2003-08-23 16:15:35 paklein Exp $ */
#ifndef _PENALTY_CONTACT_DRAG_2D_T_H_
#define _PENALTY_CONTACT_DRAG_2D_T_H_

/* base classes */
#include "PenaltyContact2DT.h"

/* direct members */
#include "InverseMapT.h"

namespace Tahoe {

/** penalty contact formulation with constant drag force */
class PenaltyContactDrag2DT: public PenaltyContact2DT
{
public:

	/** constructor */
	PenaltyContactDrag2DT(const ElementSupportT& support, const FieldT& field);

	/** initialization after constructor */
	virtual void Initialize(void);

protected:

	/** print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/** construct the residual force vector */
	virtual void RHSDriver(void);

protected:

	/** magnitude of the drag traction */
	double fDrag;

	/** gap tolerance to enable adhesion */
	double fGapTolerance;

	/** slip tolerance to enable adhesion */
	double fSlipTolerance;

	/** \name striker information */
	/** area associated with each striker node */
	/*@{*/
	dArrayT fNodalArea;
	InverseMapT fStrikerLocNumber;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _PENALTY_CONTACT_DRAG_2D_T_H_ */
