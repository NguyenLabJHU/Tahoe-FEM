//DEVELOPMENT

#ifndef _FINESCALE_T_H_ 
#define _FINESCALE_T_H_ 

#include "BCJ_MatlT.h"


namespace Tahoe {

class FineScaleT
{
public:

	enum Eqn_TypeT 		{ kVMS_BCJ, kPOWER_LAW, kPHEN };

	FineScaleT ( void ) { }
	virtual ~FineScaleT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct ( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, int =kBackward_Euler) =0;
	virtual void Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb,double delta_t=0.0	)	=0; 
  virtual void Form_RHS_F_int	(	dArrayT &F_int,double delta_t=0.0	) =0; 
	

};
} // namespace Tahoe 
#endif /* _FINESCALE_T_H_ */

