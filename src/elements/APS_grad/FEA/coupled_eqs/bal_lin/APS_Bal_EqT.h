//DEVELOPMENT

#ifndef _APS_BALEQ_T_H_ 
#define _APS_BALEQ_T_H_ 

#include "BalLinMomT.h"

namespace Tahoe {

/** APS_Bal_EqT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a coupled approach to implementation of 
 *  the Balance of Linear Momentum in weak form (i.e. the Virtual Work 
 *  Equation. **/

class APS_Bal_EqT	: public BalLinMomT
{

	public:

  	enum B_T { 
								kB, 
						   		kBgamma,
	             				kNUM_B_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

	enum VB_T {
								kN,
								knuB,
								knuNgam,
								kVB_Temp1,
								kVB_Temp2,
	             				kNUM_VB_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
								
	enum C_T { 
								kMu,
								kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

		//--------------------------------------------------------------
		
		APS_Bal_EqT 	( void ) { } 

		APS_Bal_EqT 	( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_VariableT &np1, APS_VariableT &n, 
								int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler);

		void 	Construct 		( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_VariableT &np1, APS_VariableT &n, 
								int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler); 

  		void 	Form_LHS_Keps_Kd	(	dMatrixT &Keps, dMatrixT &Kd ); // add delta_t for dynamics
  		void 	Form_RHS_F_int		(	dArrayT  &F_int ); 
		void 	Form_B_List 		( void );  // Strain Displacement Matricies
		void 	Form_VB_List 		( void );  // Strain Matricies
 		void 	Form_C_List 		( APS_MaterialT *Shear_Matl );  // Constant List

		void  	Get ( StringT &Name, FEA_dMatrixT &tensor );
		void  	Get ( StringT &Name, FEA_dScalarT &scalar );
		
		//TEMP - not needed?
		//void 	Get ( int scalar_code, FEA_dScalarT &scalar  ) { scalar = S[scalar_code]; } 

	protected:

  		FEA_dMatrix_ArrayT B; 
  		FEA_dVector_ArrayT VB; 
  		dArrayT 			C;

	protected:

		FEA_IntegrationT 		Integral;
		FEA_Data_ProcessorT 	Data_Pro; 

		double delta_t;
		int time_step;

		int n_ip, n_rows, n_cols, n_sd, n_en, n_sd_x_n_sd, n_sd_x_n_en, Time_Integration_Scheme;
  
};


} // namespace Tahoe 
#endif /* _APS_BAL_EQT_H_ */

