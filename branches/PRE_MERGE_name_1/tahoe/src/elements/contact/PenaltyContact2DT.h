/* $Id: PenaltyContact2DT.h,v 1.2 2002-06-08 20:20:19 paklein Exp $ */
/* created: paklein (12/11/1997) */

#ifndef _PENALTY_CONTACT2D_T_H_
#define _PENALTY_CONTACT2D_T_H_

/* base classes */
#include "Contact2DT.h"

class PenaltyContact2DT: public Contact2DT
{
public:

	/* constructor */
	PenaltyContact2DT(const ElementSupportT& support, const FieldT& field);

	/* writing output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the effective mass matrix */
	virtual void LHSDriver(void);

	/* construct the residual force vector */
	virtual void RHSDriver(void);
		
protected:

	double fK; // penalty "stiffness"

	/* element coords and displacements */
	dArray2DT fElCoord;
	dArray2DT fElDisp;

private:
	
	/* tracking */
	int    fnum_contact;
	double fh_max;
};

#endif /* _PENALTY_CONTACT2D_T_H_ */
