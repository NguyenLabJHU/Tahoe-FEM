/* $Id: MultiplierContactElement2DT.h,v 1.1 2002-03-18 19:24:23 rjones Exp $ */

#ifndef _MULTIPLIER_CONTACT_ELEMENT_2D_T_H_
#define _MULTIPLIER_CONTACT_ELEMENT_2D_T_H_

/* base classes */
#include "ContactElementT.h"

class MultiplierContactElement2DT: public ContactElementT
{
  public:

	/* constructor */
	MultiplierContactElement2DT
		(FEManagerT& fe_manager, XDOF_ManagerT* xdof_nodes);

	/* writing output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

	enum EnforcementParametersT { 
                                kConsistentTangent = 0 ,
                                kPenalty ,
								kGScale,
								kPScale,
								kTolP,
                                kNumEnfParameters};

	enum StatusT {	kNoP = -1,	
					kPZero,
					kPJump,
					kGapZero};

  protected:
	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* set contact status*/
	virtual void SetStatus(void);
	
	/* construct the residual force vector, called before LHS */
	virtual void RHSDriver(void);
	
	/* construct the effective mass matrix */
	virtual void LHSDriver(void);

};

#endif /* _MULTIPLIER_CONTACT_ELEMENT_2D_T_H_ */
