/* $Id: MultiplierContactElement2DT.h,v 1.10 2003-09-22 21:20:49 paklein Exp $ */
// created by : rjones 2001
#ifndef _MULTIPLIER_CONTACT_ELEMENT_2D_T_H_
#define _MULTIPLIER_CONTACT_ELEMENT_2D_T_H_

/* base classes */
#include "ContactElementT.h"

namespace Tahoe {

class MultiplierContactElement2DT: public ContactElementT
{
  public:

	/* constructor */
	MultiplierContactElement2DT(const ElementSupportT& support, const FieldT& field);

	enum EnforcementParametersT { 
                                kConsistentTangent = 0 ,
                                kPenalty ,
								kGScale,
								kPScale,
								kTolP,
								kMaterialType,
                                kNumEnfParameters};

	enum StatusT {	kNoP = -1,	
					kPZero,
					kPJump,
					kGapZero};

  protected:
	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
		 	
	/* set contact status*/
	virtual void SetContactStatus(void);
	
	/* construct the residual force vector, called before LHS */
	virtual void RHSDriver(void);
	
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

};

} // namespace Tahoe 
#endif /* _MULTIPLIER_CONTACT_ELEMENT_2D_T_H_ */
