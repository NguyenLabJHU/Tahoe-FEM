/* $Id: PenaltyContactDrag2DT.h,v 1.1 2003-08-14 05:50:41 paklein Exp $ */
#ifndef _PENALTY_CONTACT_DRAG_2D_T_H_
#define _PENALTY_CONTACT_DRAG_2D_T_H_

/* base classes */
#include "PenaltyContact2DT.h"

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

private:

	/** compute the nodal area associated with each striker node */
	void ComputeNodalArea(const ArrayT<StringT>& striker_blocks, dArrayT& nodal_area);

protected:

	/** magnitude of the drag traction */
	double fDrag;

	/** gap tolerance to enable adhesion */
	double fGapTolerance;

	/** slip tolerance to enable adhesion */
	double fSlipTolerance;

	/** area associated with each striker node */
	dArrayT fNodalArea;
};

} /* namespace Tahoe */

#endif /* _PENALTY_CONTACT_DRAG_2D_T_H_ */
