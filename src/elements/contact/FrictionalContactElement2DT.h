/* $Id: FrictionalContactElement2DT.h,v 1.2 2004-07-15 08:28:08 paklein Exp $ */
// created by : rjones 2003
#ifndef _FRICTIONAL_CONTACT_ELEMENT_2D_T_H_
#define _FRICTIONAL_CONTACT_ELEMENT_2D_T_H_

/* base classes */
#include "ContactElementT.h"

namespace Tahoe {

class FrictionalContactElement2DT: public ContactElementT
{
  public:

	/* constructor */
	FrictionalContactElement2DT
			(const ElementSupportT& support, const FieldT& field);

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
//	virtual void PrintControlData(ostream& out) const;
		 	
	/* set contact status*/
	virtual void SetContactStatus(void);
	
	/* construct the residual force vector, called before LHS */
	virtual void RHSDriver(void);
	
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT sys_type);

};

} // namespace Tahoe 
#endif /* _FRICTIONAL_CONTACT_ELEMENT_2D_T_H_ */
