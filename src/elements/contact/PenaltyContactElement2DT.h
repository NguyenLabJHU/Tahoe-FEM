/* $Id: PenaltyContactElement2DT.h,v 1.2 2002-01-28 18:43:16 dzeigle Exp $ */

#ifndef _PENALTY_CONTACT_ELEMENT_2D_T_H_
#define _PENALTY_CONTACT_ELEMENT_2D_T_H_

/* base classes */
#include "ContactElementT.h"

class PenaltyContactElement2DT: public ContactElementT
{
  public:

	/* constructor */
	PenaltyContactElement2DT(FEManagerT& fe_manager);

	/* writing output */
	virtual void WriteOutput(IOBaseT::OutputModeT mode);

        enum EnforcementParametersT { kPass = 0,
                                kPenalty ,
                                kConsistentTangent ,
                                kStatMean,
                                kStatStandDev,
                               	kAspDens, 
				kNumEnfParameters};

	 	
  protected:

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the residual force vector, called before LHS */
	virtual void RHSDriver(void);
	
	/* construct the effective mass matrix */
	virtual void LHSDriver(void);

};

#endif /* _PENALTY_CONTACT_ELEMENT_2D_T_H_ */

