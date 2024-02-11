// $Id: APS_kappa_gpT.h,v 1.2 2005-05-03 15:54:47 raregue Exp $
#ifndef _APS_kappa_gp_T_H_ 
#define _APS_kappa_gp_T_H_ 

#include "PlastT.h"

namespace Tahoe {

/** APS_kappa_gpT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for implementation of an anti-plane shear strain gradient model,
 *  with same kappa for each slip system alpha
 **/


class APS_kappa_gpT : public PlastT
{

	public:

  	enum B_d_T { 
								kB, 
						   		kBmvgam1d,
						   		kBmvgam2d,
						   		kBmvgam3d,
	             				kNUM_B_d_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
  	enum B_gradu_T { 
						   		kgrad_u,
	             				kNUM_B_gradu_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

  	enum B_gradgammap_T { 
						   		kgrad_gammap,
	             				kNUM_B_gradgammap_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)	   	             					             				
	             				
  	enum B_eps_T { 
						   		kBgamma,
						   		kBmvgam1eps,
						   		kBmvgam2eps,
						   		kBmvgam3eps,
	             				kNUM_B_eps_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)	             				
	             				
	enum VB_d_T { // Vectors 
								kVdelgam1d,
								kVdelgam2d,
								kVdelgam3d,
								kVB_d_Temp1,
								kVB_d_Temp2,
								kVB_d_Temp3,
								kVB_d_Temp4,
	             				kNUM_VB_d_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
	             				
	enum VB_eps_T { // Vectors 
								kNgam_xdy,
								kNgam_ydx,
								kVdelgam1eps,
								kVdelgam2eps,
								kVdelgam3eps,
								kVB_eps_Temp1,
								kVB_eps_Temp2,
								kVB_eps_Temp3,
								kVB_eps_Temp4,
	             				kNUM_VB_eps_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)							             				
        	
	             	
    enum V_T { 	// Vectors 
								kgammap,
								kgammap_n,
								kdel_gammap,
								km1_bar,
								km2_bar,
								km3_bar,
								kV_Temp1,
								kV_Temp2,
								kV_Temp3,
								kV_Temp4,
	             				kNUM_V_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
	enum V_state_T { // Vectors 
								kstate,
								kstate_n,
								kV_state_Temp1,
	             				kNUM_V_state_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

	enum V_out_T { 	// Vectors 
								kstressstate,
								kV_out_Temp1,
	             				kNUM_V_out_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
					
    enum S_T { 	// Scalar values
								kS_1,
								kS_2,
								kS_3,
								kgammap_curl,
								kxi_1,
								kxi_2,
								kxi_3,
								kabs_xi_1,
								kabs_xi_2,
								kabs_xi_3,
								ksignxi_1,
								ksignxi_2,
								ksignxi_3,
								kxi_curl_term,
								kzeta_1,
								kzeta_2,
								kzeta_3,
								kIV_kappa1,
								kIV_kappa1_n,
								kIV_kappa2,
								kIV_kappa2_n,
								kIV_kappa3,
								kIV_kappa3_n,
								kmag_del_gammap,
								kmag_m1,
								kmag_m2,
								kmag_m3,
								kdel_gamma1,
								kdel_gamma2,
								kdel_gamma3,
								kdel_gamma1_n,
								kdel_gamma2_n,
								kdel_gamma3_n,
								ksign_del_gamma1,
								ksign_del_gamma2,
								ksign_del_gamma3,
								kdeldel_gamma1,
								kdeldel_gamma2,
								kdeldel_gamma3,
								kA1,
								kA2,
								kA3,
								kB1,
								kB2,
								kB3,
								kS_Temp1,
								kS_Temp2,
								kS_Temp3,
								kS_Temp4,
								kS_Temp5,
								kS_Temp6,
								kS_Temp7,
								kS_Temp8,
								kS_Temp9,
								kS_Temp10,
								kS_Temp11,
								kS_Temp12,
								kS_Temp13,
								kS_Temp14,
	             				kNUM_S_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	
	enum C_T { 
								kMu,
								kgamma0_dot_1,
								kgamma0_dot_2,
								kgamma0_dot_3,
								km_rate,
								km1_x,
								km1_y,
								km2_x,
								km2_y,
								km3_x,
								km3_y,
								kl,
								kH,	
								kkappa0_1,
								kkappa0_2,
								kkappa0_3,								
								k1by3,
								kRoot3by3,
								k1byRoot3,
								ksmall,
								kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

	//--------------------------------------------------------------
	
 	APS_kappa_gpT	( void ) { }
								
	APS_kappa_gpT	( FEA_ShapeFunctionT&, FEA_ShapeFunctionT&, APS_MaterialT*, APS_VariableT&, APS_VariableT&, 
				int &fTime_Step, double fdelta_t = 0.0, int IntegrationScheme = FEA::kBackward_Euler);

	//~APS_kappa_gpT	( void ) { }

	void  	Initialize	( int &in_ip, int &in_sd, int &in_en_displ, int &in_en_plast, int &in_state, int &in_str, 
							int Initial_Time_Step=1 );

	void 	Construct 	( FEA_ShapeFunctionT&, FEA_ShapeFunctionT&, APS_MaterialT*, APS_VariableT&, APS_VariableT&, 
						int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler); 
  	void 	Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd ); 
  	void 	Form_RHS_F_int		( dArrayT &F_int ); 
	void 	Form_B_List 		( void );  // Strain Displacement Matricies
	void 	Form_VB_List 		( void );  // Strain Matricies
	void 	Form_V_S_Lists 		( APS_VariableT &np1, APS_VariableT &n); // 
 	void 	Form_C_List 		( APS_MaterialT *APS_Matl );  // Constant List

	void  	Get ( StringT &Name, FEA_dScalarT &scalar );
	void  	Get ( StringT &Name, FEA_dVectorT &vector );
	void  	Get ( StringT &Name, FEA_dMatrixT &tensor );
	//void 	Get ( int scalar_code, FEA_dScalarT &scalar  ) { scalar = S[scalar_code]; } 

	protected:

		// check these dimensions
	  	FEA_dMatrix_ArrayT    B_d, B_eps, B_gradu, B_gradgammap; 
	  	FEA_dScalar_ArrayT    S;
	  	FEA_dVector_ArrayT    V, V_state, VB_d, VB_eps, V_out; 
	  	dArrayT 			  C;
			
	protected:

		FEA_IntegrationT 	Integral;
		APS_FEA_Data_Processor_DisplT Data_Pro_Displ; 
		APS_FEA_Data_Processor_PlastT Data_Pro_Plast; 

		int n_ip, n_rows, n_cols, n_sd, n_en_displ, n_en_plast, n_sd_x_n_sd, n_sd_x_n_en_plast, Time_Integration_Scheme;
		int time_step, n_state, n_str;

		double delta_t;
};

} // namespace Tahoe 
#endif /* _APS_kappa_gp_T_H_ */

