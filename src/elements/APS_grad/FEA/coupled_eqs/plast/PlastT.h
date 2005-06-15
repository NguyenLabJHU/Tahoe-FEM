// $Id: PlastT.h,v 1.7 2005-06-15 00:06:19 raregue Exp $
#ifndef _PLAST_T_H_ 
#define _PLAST_T_H_ 

#include "StringT.h"
#include "APS_MatlT.h"
#include "APS_FEA.h"
#include "APS_EnumT.h"
#include "APS_VariableT.h"

namespace Tahoe {

/* forward declaration */
class FEA_ShapeFunctionT;

class PlastT
{
public:

	enum Eqn_TypeT 	{ kAPS_kappa_alpha, kAPS_kappa_gp, kAPS_kappa_alpha_mac };

	PlastT ( void ) { };
	virtual ~PlastT ( void ) { };

	/** Pure virtual functions */

	virtual void Construct ( FEA_ShapeFunctionT&, FEA_ShapeFunctionT&, APS_MaterialT*, APS_VariableT&, APS_VariableT&, 
						int	&fTime_Step, double  fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler ) = 0; 
	virtual void Form_LHS_Keps_Kd	(	dMatrixT &Keps, dMatrixT &Kd	)	= 0; 
	virtual void Form_RHS_F_int	(	dArrayT &F_int	) = 0;
	virtual void Initialize	(	int &in_ip, int &in_sd, int &in_en_displ, int &in_en_plast, int &in_state, int &in_str, 
								int Initial_Time_Step ) = 0;

	virtual void Get ( StringT &Name, FEA_dMatrixT &tensor ) =0;
	virtual void Get ( StringT &Name, FEA_dVectorT &vector ) =0;
	virtual void Get ( StringT &Name, FEA_dScalarT &scalar ) =0;
  	//virtual void Get ( int scalar_code, FEA_dScalarT &scalar ) =0; 

};
} // namespace Tahoe 
#endif /* _PLAST_T_H_ */

