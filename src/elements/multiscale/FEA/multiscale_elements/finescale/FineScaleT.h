// $Id: FineScaleT.h,v 1.12 2003-09-16 16:41:26 raregue Exp $
#ifndef _FINESCALE_T_H_ 
#define _FINESCALE_T_H_ 

#include "BCJ_MatlT.h"
#include "StringT.h"

namespace Tahoe {

class FineScaleT
{
public:

	enum Eqn_TypeT 		{ kVMS_BCJ, kVMS_BCJ_X, kVMS_BCJ_XX, kVMS_BCJ_RR, kVMS_EZ, kVMS_EZ2, kVMS_EZ3, kVMS_EZ4, 
											kVMS_EZ5, kPOWER_LAW, kPHEN };

	FineScaleT ( void ) { }
	virtual ~FineScaleT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct ( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
						int	&fTime_Step, double  fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler ) = 0;
	virtual void Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb	)	= 0; 
	virtual void Form_RHS_F_int	(	dArrayT &F_int	) = 0;
	virtual void Initialize	(	int &in_ip, int &in_sd, int &in_en, int Initial_Time_Step ) = 0;

	virtual void Get ( StringT &Name, FEA_dMatrixT &tensor ) =0;
	virtual void Get ( StringT &Name, FEA_dScalarT &scalar ) =0;
  	virtual void Get ( int scalar_code, FEA_dScalarT &scalar ) =0; 
  	virtual void Get ( int tensor_code, FEA_dMatrixT &tensor,int tensor_order ) =0; 

	int  Back_Stress_Type, Iso_Hard_Type;
	bool bControl_Eb, bDiagnosis_Variables, bDel_Curl_sE;

};
} // namespace Tahoe 
#endif /* _FINESCALE_T_H_ */

