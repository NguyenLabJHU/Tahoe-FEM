/* $Id: PenaltyContactElement2DT.h,v 1.22 2004-07-15 08:28:08 paklein Exp $ */
// created by : rjones 2001
#ifndef _PENALTY_CONTACT_ELEMENT_2D_T_H_
#define _PENALTY_CONTACT_ELEMENT_2D_T_H_

/* base classes */
#include "ContactElementT.h"

#include "pArrayT.h"

namespace Tahoe {

class C1FunctionT;

class PenaltyContactElement2DT: public ContactElementT
{
  public:

	/* constructor */
	PenaltyContactElement2DT(const ElementSupportT& support, const FieldT& field);
	
	/* initialize */
	virtual void Initialize(void);

	/* writing output */
	virtual void WriteOutput(void);

    enum EnforcementParametersT { 
                                kConsistentTangent = 0 ,
                                kPenalty ,
                                kMaterialType ,
								kNumEnfParameters}; 

  protected:

	/* print element group data */
//	virtual void PrintControlData(ostream& out) const;
		 	
	/* construct the residual force vector, called before LHS */
	virtual void RHSDriver(void);
	
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);
	
	/* total _real_ area of contact for each surface */
	dArrayT fRealArea; 
	dArrayT fContactArea; 
	dArrayT fPlasticArea; 

    /* penalty models */
	pArrayT<C1FunctionT*> fPenaltyFunctions;


};

} // namespace Tahoe

#endif /* _PENALTY_CONTACT_ELEMENT_2D_T_H_ */

