// $Id: PlastT.h,v 1.1 2003-07-10 17:23:53 raregue Exp $
#ifndef _PLAST_T_H_ 
#define _PLAST_T_H_ 

#include "APS_MatlT.h"
#include "StringT.h"

namespace Tahoe {

class PlastT
{
public:

	enum Eqn_TypeT 	{ kAPS_BCJ };

	PlastT ( void ) { }
	virtual ~PlastT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct ( FEA_ShapeFunctionT&, APS_MaterialT*, APS_VariableT&, APS_VariableT&, 
						int	&fTime_Step, double  fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler ) = 0; 
	virtual void Form_LHS_Keps_Kd	(	dMatrixT &Keps, dMatrixT &Kd	)	= 0; 
	virtual void Form_RHS_F_int	(	dArrayT &F_int	) = 0;
	virtual void Initialize	(	int &in_ip, int &in_sd, int &in_en, int Initial_Time_Step ) = 0;

	virtual void Get ( StringT &Name, FEA_dMatrixT &tensor ) =0;
	virtual void Get ( StringT &Name, FEA_dScalarT &scalar ) =0;
  	virtual void Get ( int scalar_code, FEA_dScalarT &scalar ) =0; 

};
} // namespace Tahoe 
#endif /* _PLAST_T_H_ */

