
//DEVELOPMENT

#ifndef _V_BALLINMOM_T_H_ 
#define _V_BALLINMOM_T_H_ 

#include "StringT.h"
#include "Shear_MatlT.h"
#include "APS_MatlT.h"
#include "APS_V_FEA.h"
#include "APS_EnumT.h"
#include "APS_VariableT.h"

namespace Tahoe {

class V_BalLinMomT
{
public:

	enum Eqn_TypeT { kAPS_V_Bal_Eq };

	V_BalLinMomT ( void ) { }
	virtual ~V_BalLinMomT ( void ) { }

	/** Pure virtual functions */

	virtual void Construct ( int& nipsurf, int& nensurf, FEA_ShapeFunctionT&, FEA_ShapeFunctionT&, APS_MaterialT*, 
							APS_MaterialT*, 
							APS_VariableT&, APS_VariableT&, 
							int	&fTime_Step, double fdelta_t = 0.0, int =FEA::kBackward_Euler) =0;
	virtual void Form_LHS_Keps_Kd	(	dMatrixT &Keps, dMatrixT &Kd	)	=0; 
  	virtual void Form_RHS_F_int	( dArrayT  &F_int, APS_VariableT &npt ) =0; 
  	virtual void Form_LHS_Kd_Surf ( dMatrixT &Kd_face, FEA_SurfShapeFunctionT &SurfShapes ) =0;
  	virtual void Form_RHS_F_int_Surf ( dArrayT  &F_int_face, APS_VariableT &npt, double &wght  ) =0; 
	virtual void Get ( StringT &Name, FEA_dScalarT &scalar ) =0;
	virtual void Get ( StringT &Name, FEA_dVectorT &vector ) =0;
	virtual void Get ( StringT &Name, FEA_dMatrixT &tensor ) =0;
  	//virtual void Get ( int scalar_code, FEA_dScalarT &scalar ) =0; 

};
} // namespace Tahoe 
#endif /* _V_BALLINMOM_T_H_ */

