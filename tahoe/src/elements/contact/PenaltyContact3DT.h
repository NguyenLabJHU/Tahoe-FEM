/* $Id: PenaltyContact3DT.h,v 1.5 2002-11-30 16:41:27 paklein Exp $ */
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

	/* writing output */
	virtual void WriteOutput(void);
	 	
protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/* construct the residual force vector */
	virtual void RHSDriver(void);
	
private:

	/* set surface normal derivative matrix */
	void Set_dn_du(const dArray2DT& curr_coords, dMatrixT& dn_du) const;
	
protected:

	double fK; // penalty "stiffness"

	/* element coords and displacements */
	dArray2DT fElCoord;
	dArray2DT fElDisp;
	
	/* work space */
	dMatrixT fdc_du;
	dMatrixT fdn_du;
	dMatrixT fM1;
	dMatrixT fM2;
	dArrayT  fV1;
	
private:

	AutoArrayT<double> fDists; // set during RHS, used for LHS

	/* tracking */
	int    fnum_contact;
	double fh_max;
};

} // namespace Tahoe

#endif /* _PENALTY_CONTACT3D_T_H_ */
