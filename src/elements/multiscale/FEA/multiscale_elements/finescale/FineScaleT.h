// $Id: FineScaleT.h,v 1.8 2003-02-03 04:40:27 paklein Exp $
#ifndef _FINESCALE_T_H_ 
#define _FINESCALE_T_H_ 

#include "BCJ_MatlT.h"

namespace Tahoe {

class FineScaleT
{
public:

	enum Eqn_TypeT 		{ kVMS_BCJ, kVMS_EZ, kVMS_EZ2, kVMS_EZ3, kVMS_EZ4, 
											kVMS_EZ5, kPOWER_LAW, kPHEN };

	FineScaleT ( void ) { }
	virtual ~FineScaleT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct ( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
								double  fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler ) = 0;
	virtual void Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb	)	= 0; 
	virtual void Form_RHS_F_int	(	dArrayT &F_int	) = 0; 

};
} // namespace Tahoe 
#endif /* _FINESCALE_T_H_ */

