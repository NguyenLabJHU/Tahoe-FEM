// $Id: VMS_BCJT.h,v 1.11 2003-04-23 23:34:25 creigh Exp $
#ifndef _VMS_BCJT_H_ 
#define _VMS_BCJT_H_ 

#include "FineScaleT.h"

namespace Tahoe {

/** VMS_BCJTT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a Variational Multi-Scale (VMS) approach to implementation of Sandia's
 *  BCJT Model.  A dual field formulation u^alpha and u^beta is used. See
 *  Creighton et. al for more. 
 *  Collaboration:  Sandia National Laboratory and the University of Michigan **/


class VMS_BCJT : public FineScaleT
{

	public:

  	enum B_T { 
								kB, 
						   	kB_1hat,   
							 	kB05_tau_3hat,
							 	kB06_3hat,
								
								kBa_Cb_3hat,
							 	kBb_Cb_3hat,
								kBa_Cb_tau_3hat,
							 	kBb_Cb_tau_3hat,
								kBa_Cbi_3hat,
								kBb_Cbi_3hat,
								kBa_Cbi_tau_3hat,
								kBb_Cbi_tau_3hat,
								kBa_H,
								kBb_H,
								kBa_S,
								kBb_S,
								kBa_Z,
								kBb_Z,
								kBa_Zeta_Jb,
								kBb_Zeta_Jb,
								kBa_Zeta_Fbi,
								kBb_Zeta_Fbi,
								kBa_Zeta_Fbi_reg,
								kBb_Zeta_Fbi_reg,
								kBa_Zeta_Fbi_tau,
								kBb_Zeta_Fbi_tau,
								kBa_Zeta_curl_sE_T,
								kBb_Zeta_curl_sE_T,
								kBa_DEV_H, 
								kBb_DEV_H, 

								kBa_Kappa_tau_3hat,   // Iso Hardening
								kBa_Kappa_3hat,
								kBa_Kappa,

								kB_Temp0,
								kB_Temp1,
	             	kNUM_B_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

								
    enum A_T { 
						   	kgrad_ub, 
						   	kF,   		// for Iso hard
						   	kFT,  		// for Iso hard
						   	kFa,  		// for Iso hard
						   	kFaT, 		// for Iso hard
							 	kFa_n,		// for Iso hard
						   	kFai,
						   	kFa_dot, 	// for expr G2.
						   	kFb,			//
								kFbT,
							 	kCa_n,
							 	kLa,			// Diagnostics
							 	kDa,			// Diagnostics
							 	kla,			// Diagnostics
							 	kda,			// Diagnostics
							 	kNda,			// Diagnostics
							 	kDa_m,
							 	kDa_mp1,
								kF_sharp,
								kF_sharp_T,
								kN,
								kNea,		// Direction of Ea
						   	kNeb,		// Direction of Eb
								kA1,
								kA2,
								kA3,
								kA1T,
								kA2T,
								kA3T,		
						   	kC,			// Diagnostics
						   	kCa,			
						   	kCb,			
						   	kCbi,			
						   	kE,			// Diagnostics
						   	kEa,			
						   	kEa_n,		
						   	kEa_dot,	
						   	kEb,			
						   	kSigma, // Diagnostics	
						   	kS,				
						   	kS_tilde,			// Diagnostics
						   	kS_hat,				// Diagnostics
						   	kFbi,					// for Back Stress
						   	kH_bar,				// for Back Stress
						   	kDEV_H,				// for Back Stress
						   	kZeta,				// for Back Stress
						   	kZetaT,				// for Back Stress
						   	kSym_Zeta,		// for Back Stress
						   	kcurl_sE,			// for Back Stress
								kA4,					// for Back Stress
								kA4T,					// for Back Stress
						   	kMu_c_l_Jb_curl_sE_T,	// for Back Stress
						   	kR2T,					// for Back Stress
						   	kFbiT,				// for Back Stress
								kG2,
							 	kFb_n,				// for CCba GammaY 
							 	kCb_n,				// for CCba GammaY 
							 	kEb_n,				// for CCba GammaY 
							 	kEb_dot,			// for CCba GammaY 
							 	kN_Eb_dot,		// for CCba GammaY 
								kA_Temp0,
	             	kNUM_A_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

								
    enum T4_T { 	//	4th order tensors in reduced form 

								kCC,					//	Doubled upper case represent black bold font 
								kCC1,					//	Used for CCvmst 
								kCC2,					//	Used for CCvmst 
								kCC_ba,				//	Used for CCvmst 
								kCC_gamma,		//	Used for CCvmst 
								kCC_psi,			//	Used for CCvmst 
								kEb_o_Eb,			//	Used for CCvmst 
								kEEn,					//	Used for CCvmst 

								kMM,
								kPP,
								kMag_DEV_H_PP,
								kN_o_N,		
								kN_o_Nea,			// Used for kappa
								kCbi_o_H,			// Outer product of (Cb)^-1 and S  
								kCbi_o_Cb,
								kAAe_2bar,			// in Linearization of Zeta
								kAAf_2bar,			// in Linearization of Zeta
								kEE_2bar,				// in Linearization of Zeta
								kZ_o_F_sharp_T,	// in Linearization of Zeta
								kZ_o_1,					// in Linearization of Zeta
								kT4_Temp0,
	             	kNUM_T4_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

								
    enum S_T { 	// Scalar values

								kMag_DEV_H,
								kMag_Ea,
						   	kMag_Ea_dot,
								kBeta,
								kBeta2,
								kBetaK,
								kKappa_bar,
								kMacaulay_Sinh_Beta,
								kCosh_Beta, // No point in Macaulay Cosh(x) min point is 1.0
								kCb_i_H,		// Inner product of C^beta and H 
								kIV_Alpha,	// Internal Variable called Alpha
								kIV_Alpha_n,
								kJb,
								kJ,
								kpLamda, // for expr G2.Plastic Lamda as in:  S = pLamda*Da (Visco)

								kRho_Mag_Eb,
								kMag_Eb,
								kPi_Tanh,
								kPsi,
								kPsi1,
								kEta1,
								kEta2,
								kEta1_plus_Eta2,
								kGammaY_bar_Eta1_plus_Eta2,
								kLamda_ba,
								kMu_ba,
								kMag_Eb_dot,	
								kGammaY,			
								kGammaY_bar,			

								kMag_Da,			// Diagnostics
								kArcBeta,			// Diagnostics
								kArcSinh,			// Diagnostics
								kS_Temp0,
								kS_Temp1,
	             	kNUM_S_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

								
   enum C_T { 
								kLamda,
								kMu,
								kf,
								kV,
								kY,
								kc,
								kl,
								kMu_c_l,
								kKappa,
								kH,	

								kLamda1,
								kMu1,
								kLamda2,
								kMu2,
								kPi,
								kRho,
								kGamma_b,
								k2Gamma_b,
								kAlphaY,
								kAlphaY_plus_1,

								kNeg_dt_Root3by2_f,
								kRoot3by2byV,
								k1by2,
								k1by3,
								kRoot3by2,
								kRoot2by3,
								k2by3,
	             	kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

		enum NP_T {
								ksE_mat,
								kFbi_hat_mat,
								kIdentity,
								kNUM_NP_TERMS };


		enum BCK_T {
								kNo_Back_Stress,
								kSteinmann,
								kExp_Back };

		enum ISO_T {
								kNo_Iso_Hard,
								kMethod1,
								kMethod2,
								kMethod3,
								kExp_Iso };

	//--------------------------------------------------------------
	
 	VMS_BCJT	(	) { }
								
	VMS_BCJT	( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
							int &fTime_Step, double fdelta_t = 0.0, int IntegrationScheme = FEA::kBackward_Euler);

	void  Initialize	( int &in_ip, int &in_sd, int &in_en, int Initial_Time_Step=1 );

	void 	Construct 	( FEA_ShapeFunctionT&, VMF_MaterialT*, VMS_VariableT&, VMS_VariableT&, 
								 			int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler); 

  void 	Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb ); 
  void 	Form_RHS_F_int	(	dArrayT &F_int ); 
	void 	Form_B_List 		( void );  // Strain Displacement Matricies
	void 	Form_A_S_Lists 	( VMS_VariableT &np1, VMS_VariableT &n); // BCDE ---> A
	void 	Form_T4_List 		( void );  // Tensors generated by order reduction of 4th order Tensor
 	void 	Form_C_List 		( VMF_MaterialT *BCJ_Matl );  // Constant List

	void  Get ( StringT &Name, FEA_dScalarT &scalar );
	void  Get ( StringT &Name, FEA_dMatrixT &tensor );
	void 	Get ( int scalar_code, FEA_dScalarT &scalar  ) { scalar = S[scalar_code]; } 
	void 	Get ( int tensor_code, FEA_dMatrixT &tensor, int tensor_order ) 		
							{ tensor = (tensor_order==2) ? A[tensor_code] : T4[tensor_code]; } 

	void 	CCba_Scalars		( void );	
	void	Form_CC					( void ); // Fills CC with either classic CC or CCba
	
	void 	Get_Iso_Hard_Kappa_bar ( void );
	void 	Get_Back_Stress	( void );
	void  Form_del_Zeta_B	( void );  // Calculates B[kBa_Z],  B[kBb_Z]	



	protected:

  	FEA_dMatrix_ArrayT    B; 	// Dimension: n_sd x n_sd*n_en , 		e.g. 9x24
		FEA_dMatrix_ArrayT    A; 	// Dimension: n_sd x n_sd , 					e.g. 3x3 
  	FEA_dMatrix_ArrayT   T4; 	// Dimension: n_sd*n_sd x n_sd*n_sd, e.g. 9x9 Matricies 
															// (Tensors generated from 3x3 outer 3x3 reduced order)
  	FEA_dScalar_ArrayT    S; 
  	dArrayT 			      	C;

		// NP for Nodal Point, This special array has values at node points NOT integration points
		FEA_dMatrix_ArrayT    NP; 	
	 					
	protected:

		FEA_IntegrationT 		Integral;
		FEA_Data_ProcessorT Data_Pro; 

		int n_ip, n_rows, n_cols, n_sd, n_en, n_sd_x_n_sd, n_sd_x_n_en, Time_Integration_Scheme;
		int time_step;

		double delta_t;
};

} // namespace Tahoe 
#endif /* _VMS_BCJT_H_ */

