
//DEVELOPMENT

#ifndef _MFGP_BALLINMOM_T_H_ 
#define _MFGP_BALLINMOM_T_H_ 

#include "StringT.h"
#include "GRAD_MR_Plast_MatT.h" //include the Gradient plas model
#include "MFGP_MFA.h"
#include "MFGP_EnumT.h"
#include "MFGP_VariableT.h"

namespace Tahoe 
{

class MFGP_BalLinMomT
{
public:

	enum Eqn_TypeT { kMFGP_Bal_Eq };

	MFGP_BalLinMomT ( void ) { }
	virtual ~MFGP_BalLinMomT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct (ShapeFunctionT&, ShapeFunctionT&, GRAD_MR_Plast_MaterialT*, 
							GRAD_MR_Plast_MaterialT*, 
							MFGP_VariableT&, MFGP_VariableT&, 
							int	&fTime_Step, double fdelta_t = 0.0, int =FEA::kBackward_Euler) =0;
	virtual void Form_LHS_Klambda_Ku	(	dMatrixT &Klambda, dMatrixT &Ku	)	=0; 
  	virtual void Form_RHS_F_int	( dArrayT  &F_int, MFGP_VariableT &npt ) =0; 
	
};
} // namespace Tahoe 
#endif /* _MFGP_BALLINMOM_T_H_ */

