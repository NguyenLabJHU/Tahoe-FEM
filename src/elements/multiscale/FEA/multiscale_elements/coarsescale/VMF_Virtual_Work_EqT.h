//DEVELOPMENT

#ifndef _VMF_VWEQ_T_H_ 
#define _VMF_VWEQ_T_H_ 

#include "CoarseScaleT.h"
#include "Iso_MatlT.h"

namespace Tahoe {

/** VMF_Virtual_Work_EqT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a Variational Multi-Scale (VMS) approach to implementation of 
 *  the Conservation of Linear Momentum in weak form (i.e. the Virtual Work 
 *  Equation. **/

class VMF_Virtual_Work_EqT	: public CoarseScaleT
{

	public:

  	enum B_T { 
						   	kB_1hat,   
							 	kBI_tau_3hat,
							 	kBbII_2hat,
							 	kBaII_3hat,
							 	kBaIII_2bar,
						   	kB_Temp0,
	             	kNUM_B_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

	    enum A_T { 
						   	kF,  // NOTE: kF != VMS::kF 
						   	kFb,
						   	kgrad_ub,
						   	kFbT,
						   	kF_sharp,
						   	kCb,
						   	kEb,
						   	kNeb,
						   	kEb_tilde,
						   	kS,
						   	kSigma,
							 	kFb_n,				// for CCba GammaY 
							 	kCb_n,				// for CCba GammaY 
							 	kEb_n,				// for CCba GammaY 
							 	kEb_dot,			// for CCba GammaY 
							 	kN_Eb_dot,		// for CCba GammaY 
						   	kA_Temp0,
            		kNUM_A_TERMS };  

    	enum T4_T { 
								kCC,
								kcc_b,
								kCC1,
								kCC2,
								kCC_ba,				//	Used for CCvmst 
								kCC_gamma,		//	Used for CCvmst 
								kCC_psi,			//	Used for CCvmst 
								kEb_o_Eb,			//	Used for CCvmst 
								kEEn,					//	Used for CCvmst 
						   	kT4_Temp0,
	             	kNUM_T4_TERMS }; 

	    	enum S_T { 
								kJ,
								kJb,
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
								kMag_Eb_dot,	// CCba GammaY
								kGammaY,			// CCba GammaY
								kGammaY_bar,	// CCba GammaY
	             	kNUM_S_TERMS };  

   enum C_T { 
								kLamda,
								kMu,

								kLamda1,
								kMu1,
								kLamda2,
								kMu2,
								kGamma_b,			
								k2Gamma_b,			
								kAlphaY,			
								kAlphaY_plus_1,			
								kPi,
								kRho,

	             	kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)

		//--------------------------------------------------------------
		
		VMF_Virtual_Work_EqT 	( void ) { } 

		VMF_Virtual_Work_EqT 	( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
														int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler);

		void 	Construct 			( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
														int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler); 

  	void 	Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb ); // add delta_t for dynamics
  	void 	Form_RHS_F_int	(	dArrayT  &F_int ); 
		void 	Form_B_List 		( void );  // Strain Displacement Matricies
		void 	Form_A_S_Lists 	( VMS_VariableT &np1, VMS_VariableT &n,int Integration_Scheme=FEA::kBackward_Euler ); // BCDE ---> A 
 		void 	Form_C_List 		( VMF_MaterialT *Iso_Matl );  // Constant List

		void 	CCba_Scalars		( void );	
		void	Form_CC					( void ); // Fills CC with either classic CC or CCba

		void  Get ( StringT &Name, FEA_dMatrixT &tensor );
		void  Get ( StringT &Name, FEA_dScalarT &scalar );
		void 	Get ( int scalar_code, FEA_dScalarT &scalar  ) { scalar = S[scalar_code]; } 
		void 	Get ( int tensor_code, FEA_dMatrixT &tensor, int tensor_order ) 		
								{ tensor = (tensor_order==2) ? A[tensor_code] : T4[tensor_code]; } 

	protected:

  	FEA_dMatrix_ArrayT B; 
  	FEA_dMatrix_ArrayT A; 
  	FEA_dMatrix_ArrayT T4; 
  	FEA_dScalar_ArrayT S; 
  	dArrayT 			     C;

	protected:

		FEA_IntegrationT 		Integral;
		FEA_Data_ProcessorT Data_Pro; 

		double delta_t;
		int time_step;

		int n_ip, n_rows, n_cols, n_sd, n_en, n_sd_x_n_sd, n_sd_x_n_en, Time_Integration_Scheme;
  
};


} // namespace Tahoe 
#endif /* _VMS_BCJT_H_ */

