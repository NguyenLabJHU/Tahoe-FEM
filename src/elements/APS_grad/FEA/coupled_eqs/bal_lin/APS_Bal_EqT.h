// $Id: APS_Bal_EqT.h,v 1.10 2003-10-08 17:45:13 raregue Exp $
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

  	enum B_d_T { 
								kB,
								kB_surf, 
	             				kNUM_B_d_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

  	enum B_eps_T {  
						   		kBgamma,
	             				kNUM_B_eps_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

  	enum B_gradu_T { 
						   		kgrad_u,
						   		kgrad_u_surf,
	             				kNUM_B_gradu_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)


	enum VB_d_T {
								kN,
								knuB,
								kVB_d_Temp1,
								kVB_d_Temp2,
	             				kNUM_VB_d_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

	enum VB_eps_T {
								knuNgam,
								kVB_eps_Temp1,
								kVB_eps_Temp2,
	             				kNUM_VB_eps_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

	             					
	enum V_T {
								knueps,
								keps,
								//kgrad_u,
								kgammap,
								kV_Temp1,
								kV_Temp2,
	             				kNUM_V_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
	enum VS_T {
								kVS_Temp1,
								kVS_Temp2,
	             				kNUM_VS_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
	enum S_T {
								knuepsgradu,
								knuepseps,
								kS_Temp1,
	             				kNUM_S_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
								
	enum C_T { 
								kMu,
								km1,
								km2,
								kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

		//--------------------------------------------------------------
		
		APS_Bal_EqT 	( void ) { } 

		APS_Bal_EqT 	( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_MaterialT *APS_Matl, 
						APS_VariableT &np1, APS_VariableT &n, 
						int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler);

		void 	Construct 	( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_MaterialT *APS_Matl, 
							APS_VariableT &np1, APS_VariableT &n, 
							int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler); 

  		void 	Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd ); // add delta_t for dynamics
  		void 	Form_RHS_F_int		( dArrayT  &F_int, APS_VariableT &npt ); 
  		void 	Form_LHS_Kd_Surf	( dMatrixT &Kd, FEA_SurfShapeFunctionT &SurfShapes, iArrayT& face_equations); // add delta_t for dynamics
  		void 	Form_RHS_F_int_Surf	( dArrayT  &F_int, APS_VariableT &npt, double &wght, iArrayT& face_equations  ); 
		void 	Form_B_List 		( void );  // Strain Displacement Matricies
		void 	Form_VB_List 		( void );  // Strain Matricies
		void 	Form_V_S_List 		( APS_VariableT &npt );  // vectors
 		void 	Form_C_List 		( APS_MaterialT *Shear_Matl,  APS_MaterialT *APS_Matl );  // Constant List

		void  	Get ( StringT &Name, FEA_dMatrixT &tensor );
		void  	Get ( StringT &Name, FEA_dVectorT &vector );
		void  	Get ( StringT &Name, FEA_dScalarT &scalar );
		//void 	Get ( int scalar_code, FEA_dScalarT &scalar  ) { scalar = S[scalar_code]; } 

	protected:

  		FEA_dMatrix_ArrayT B_d, B_eps, B_gradu; 
  		FEA_dVector_ArrayT VB_d, VB_eps, V, VS;
  		FEA_dScalar_ArrayT S; 
  		dArrayT 			C;

	protected:

		FEA_IntegrationT 		Integral;
		FEA_SurfIntegrationT 	SurfIntegral;
		APS_FEA_Data_ProcessorT Data_Pro;
		APS_FEA_Data_Processor_SurfT Data_Pro_Surf; 

		double delta_t;
		int time_step;

		int n_ip, n_rows_vector, n_rows_matrix, n_cols_matrix, n_sd, n_en, n_sd_x_n_sd, 
			n_sd_x_n_en, Time_Integration_Scheme;
  
};


} // namespace Tahoe 
#endif /* _APS_BAL_EQT_H_ */

