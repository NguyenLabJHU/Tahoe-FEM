/* $Id: PenaltyContact2DT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (12/11/1997)                                          */
/* penalty based contact element. currently not well-written for          */
/* explicit dynamics. see code at end of .cpp file for what should        */
/* be added.                                                              */

#ifndef _PENALTY_CONTACT2D_T_H_
#define _PENALTY_CONTACT2D_T_H_

/* base classes */
#include "Contact2DT.h"

class PenaltyContact2DT: public Contact2DT
{
public:

	/* constructor */
	PenaltyContact2DT(FEManagerT& fe_manager);

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
