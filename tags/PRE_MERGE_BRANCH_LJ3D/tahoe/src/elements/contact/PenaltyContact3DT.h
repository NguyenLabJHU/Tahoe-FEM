/* $Id: PenaltyContact3DT.h,v 1.6 2003-03-02 18:56:12 paklein Exp $ */
/* created: paklein (02/09/2000) */
#ifndef _PENALTY_CONTACT3D_T_H_
#define _PENALTY_CONTACT3D_T_H_

/* base classes */
#include "Contact3DT.h"

namespace Tahoe {

class PenaltyContact3DT: public Contact3DT
{
public:

	/* constructor */
	PenaltyContact3DT(const ElementSupportT& support, const FieldT& field);

protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

	/* construct the residual force vector */
	virtual void RHSDriver(void);

protected:

	/** penalty "stiffness" */
	double fK;

	/** \name element coords and displacements */
	/*@{*/
	dArray2DT fElCoord;
	dArray2DT fElRefCoord;
	dArray2DT fElDisp;
	/*@}*/
	
	/** \name work space */
	/*@{*/
	dMatrixT fdc_du;
	dMatrixT fdn_du;
	dMatrixT fM1;
	dMatrixT fM2;
	dArrayT  fV1;
	/*@}*/
};

} // namespace Tahoe

#endif /* _PENALTY_CONTACT3D_T_H_ */
