/* $Id: PenaltyContact3DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (02/09/2000)                                          */
/* penalty based contact element                                          */

#ifndef _PENALTY_CONTACT3D_T_H_
#define _PENALTY_CONTACT3D_T_H_

/* base classes */
#include "Contact3DT.h"

class PenaltyContact3DT: public Contact3DT
{
public:

	/* constructor */
	PenaltyContact3DT(FEManagerT& fe_manager);

	/* writing output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);
	 	
protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the effective mass matrix */
	virtual void LHSDriver(void);

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

#endif /* _PENALTY_CONTACT3D_T_H_ */
