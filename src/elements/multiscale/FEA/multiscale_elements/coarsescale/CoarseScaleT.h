// $Id: CoarseScaleT.h,v 1.4 2003-02-03 04:40:26 paklein Exp $
#ifndef _COARSESCALE_T_H_ 
#define _COARSESCALE_T_H_ 

#include "Iso_MatlT.h"

namespace Tahoe {

class CoarseScaleT
{
public:

	enum Eqn_TypeT { kVMF_Virtual_Work_Eq, kLDV, kStraight };
				

	CoarseScaleT ( void ) { }
	virtual ~CoarseScaleT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct ( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
														int =FEA::kBackward_Euler) =0;

	virtual void Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb,double delta_t=0.0	)	=0; 
  virtual void Form_RHS_F_int	(	dArrayT &F_int, double delta_t=0.0	) =0; 
  virtual void Sigma 					(	FEA_dMatrixT &sigma ) =0; 
	

};
} // namespace Tahoe 
#endif /* _FINESCALE_T_H_ */

