// $Id: MFGP_Bal_EqT.h,v 1.1 2004-06-18 19:20:15 kyonten Exp $
#ifndef _MFGP_BAL_EQ_T_H_ 
#define _MFGP_BAL_EQ_T_H_ 

#include "MFGP_BalLinMomT.h"

namespace Tahoe 
{

/** MFGP_Bal_EqT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a coupled approach to implementation of 
 *  the Balance of Linear Momentum in weak form 
 **/

class MFGP_Bal_EqT	: public MFGP_BalLinMomT
{

	public:

  	enum B1_d_T {  //include B1 corresponding to disp
								kB1,
	             				kNUM_B1_d_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
	enum B3_d_T {  //include B3 corresponding to disp (gradient term)
								kB3,
	             				kNUM_B3_d_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	             				
  	enum phi_lam_T {  //include shape function corresponding to plast multiplier
						   		kphi,
	             				kNUM_phi_lam_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
  	
  	enum B4_lam_T {  //include B4 corresponding to plast multiplier
						   		kB4,
	             				kNUM_B4_lam_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
							
	enum C_T {  //fModuli from the constitutive model??
								kModuli,
								kNUM_C_TERMS };  // <-- Use for loops and count (KEEP THIS ONE LAST!!)
	// C matrices needed to be defined separately??
	
		//--------------------------------------------------------------
		
		MFGP_Bal_EqT 	( void ) { } 

		MFGP_Bal_EqT 	( ShapeFunctionT &Shapes_displ, ShapeFunctionT &Shapes_plast, 
						GRAD_MR_Plast_MaterialT *GRAD_MR_Plast_Mat, 
						MFGP_VariableT &np1, MFGP_VariableT &n, 
						int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler);

		void 	Construct 	(ShapeFunctionT &Shapes_displ, ShapeFunctionT &Shapes_plast, 
							MFGP_MaterialT *Shear_Matl, 
							MFGP_VariableT &np1, MFGP_VariableT &n, 
							int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme=FEA::kBackward_Euler); 

  		void 	Form_LHS_Klambda_Ku( dMatrixT &Klambda, dMatrixT &Ku); // add delta_t for dynamics
  		void 	Form_RHS_F_int( dArrayT  &F_int, MFGP_VariableT &npt ); 
		void 	Form_B_List( void );  // Strain Displacement Matricies
 		void 	Form_C_List( GRAD_MR_Plast_MaterialT *GRAD_MR_Plast_Mat);  // Constant List


	protected:
		dMatrixT    B1_d, B3_d; 
		double B4_lam, phi_lam; //dimension??
  		dMatrixT Cuu1, Cuu2, Culam1, Culam2;
  		dMatrixT Clamu1, Clamu2;
  		double Clamlam1, Clamlam2;  
  		
	protected:

		MFA_IntegrationT 		Integral;
		MFGP_MFA_Data_Processor_DisplT Data_Pro_Displ;
		MFGP_MFA_Data_Processor_PlastT Data_Pro_Plast;

		double delta_t;
		int time_step;

		int n_ip, n_rows_vector, n_rows_matrix, n_cols_matrix, n_sd, n_en_displ, n_en_plast, n_sd_x_n_sd, 
			n_sd_x_n_en_displ, n_sd_x_n_en_plast, Time_Integration_Scheme;
  
};


} // namespace Tahoe 
#endif /* _MFGP_BAL_EQT_H_ */

