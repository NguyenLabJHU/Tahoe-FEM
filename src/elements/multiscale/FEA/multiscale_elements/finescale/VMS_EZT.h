// $Id: VMS_EZT.h,v 1.5 2003-03-17 22:05:32 creigh Exp $
#ifndef _VMS_EZ_T_H_ 
#define _VMS_EZ_T_H_ 

#include "FineScaleT.h"

namespace Tahoe {

/** VMS_EZT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a Variational Multi-Scale (VMS) approach to implementation of Sandia's
 *  BCJ Model.  A dual field formulation u^alpha and u^beta is used. See
 *  Creighton et. al for more. 
 *  Collaboration:  Sandia National Laboratory and the University of Michigan **/


class VMS_EZT : public FineScaleT
{

	public:

  	enum B_T { 
						   	kB_1hat,   
	             	kNUM_B_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

    enum A_T { 
						   	kgrad_ua,  
						   	kgrad_ub,  
						   	kG2,
	             	kNUM_A_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

  	enum C_T { 
						   	kNeg_Alpha,   
	             	kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

		//--------------------------------------------------------------
	
 	VMS_EZT	(	) { }
								
	VMS_EZT	( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
						int &fTime_Step, double fdelta_t = 0.0, int IntegrationScheme = FEA::kBackward_Euler);

	void  Initialize 	( int &in_ip, int &in_sd, int &in_en, int Initial_Time_Step );
	
	void 	Construct 	( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
											int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler); 

  void 	Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb ); 
  void 	Form_RHS_F_int	(	dArrayT &F_int ); 
	void 	Form_B_List 		( void );  // Strain Displacement Matricies
	void 	Form_A_S_Lists 	( VMS_VariableT &np1, VMS_VariableT &n ); // BCDE ---> A 
 	void 	Form_C_List 		( VMF_MaterialT *BCJ_Matl );  // Constant List

	void  Get ( StringT &Name, FEA_dScalarT &scalar ) { scalar = 0.0; cout<<"EZ: Unknown scalar '"<<Name<<"' requested.\n"; } 
	void  Get ( StringT &Name, FEA_dMatrixT &tensor );
	void 	Get ( int scalar_code, FEA_dScalarT &scalar  )  { scalar_code=scalar.Length(); } // { scalar = S[scalar_code]; } 
	void 	Get ( int tensor_code, FEA_dMatrixT &tensor, int tensor_order ) 		
							{ tensor = A[tensor_code]; tensor_order++; } 
							//{ tensor = (tensor_order==2) ? A[tensor_code] : T4[tensor_code]; } 

	protected:

  	FEA_dMatrix_ArrayT    B; // Dimension: n_sd x n_sd*n_en , 		e.g. 9x24
		FEA_dMatrix_ArrayT    A; // Dimension: n_sd x n_sd , 					e.g. 3x3 
   	dArrayT 			      	C;
 
	protected:

		FEA_IntegrationT 		Integral;
		FEA_Data_ProcessorT Data_Pro; 

		int n_ip, n_rows, n_cols, n_sd, n_en, n_sd_x_n_sd, n_sd_x_n_en;
  
};

} // namespace Tahoe 
#endif /* _VMS_EZT_H_ */


