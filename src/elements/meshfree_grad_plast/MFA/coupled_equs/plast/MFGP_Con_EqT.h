// $Id: MFGP_Con_EqT.h
#ifndef _MFGP_CON_EQ_T_H_ 
#define _MFGP_CON_EQ_T_H_ 

#include "MFGP_PlastT.h"

namespace Tahoe 
{

/** MFGP_Con_EqT: This class contains methods which build stiffness matricies 
 *  and formulate the non-linear Newton-Raphson equations Kd = -R
 *  for a coupled approach to implementation of 
 *  the Consistency Condition in weak form 
 **/


class MFGP_Con_EqT : public MFGP_PlastT
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
 	MFGP_Con_EqT	( void ) { }
								
	MFGP_Con_EqT	( ShapeFunctionT&, ShapeFunctionT&, GRAD_MR_Plast_MaterialT*, MFGP_VariableT&, MFGP_VariableT&, 
				int &fTime_Step, double fdelta_t = 0.0, int IntegrationScheme = FEA::kBackward_Euler);

	void  	Initialize	( int &in_ip, int &in_sd, int &in_en_displ, int &in_en_plast, int &in_state, int &in_str, 
							int Initial_Time_Step=1 );

	void 	Construct 	( ShapeFunctionT&, ShapeFunctionT&, GRAD_MR_Plast_MaterialT*, MFGP_VariableT&, MFGP_VariableT&, 
						int &fTime_Step, double fdelta_t = 0.0, int Integration_Scheme = FEA::kBackward_Euler); 
  	void 	Form_LHS_Ku_Klambda	( dMatrixT &Ku, dMatrixT &Klambda ); 
  	void 	Form_RHS_F_int		( dArrayT &F_int ); 
	void 	Form_B_List 		( void );  // Strain Displacement Matricies
 	void 	Form_C_List 		( GRAD_MR_Plast_MaterialT *GRAD_MR_Plast_Mat);  // Constant List

	
	protected:

		// check the dimensions!!
	  	dMatrixT    B1_d, B3_d; 
		double B4_lam, phi_lam; //dimension??
	  	dMatrixT Cuu1, Cuu2, Culam1, Culam2;
  		dMatrixT Clamu1, Clamu2;
  		double Clamlam1, Clamlam2;  
			
	protected:

		MFA_IntegrationT 	Integral;
		MFGP_MFA_Data_Processor_DisplT Data_Pro_Displ; 
		MFGP_MFA_Data_Processor_PlastT Data_Pro_Plast; 

		int n_ip, n_rows, n_cols, n_sd, n_en_displ, n_en_plast, n_sd_x_n_sd, n_sd_x_n_en_plast, Time_Integration_Scheme;
		int time_step, n_state, n_str;

		double delta_t;
};

} // namespace Tahoe 
#endif /* _MFGP_CON_EQ_T_H_ */

