// $Id: MFGP_Bal_EqT.h,v 1.6 2004-08-19 18:21:55 raregue Exp $
#ifndef _MFGP_BAL_EQ_T_H_ 
#define _MFGP_BAL_EQ_T_H_ 

#include "MFGP_BalLinMomT.h"
#include "MFGP_MFA_Data_Processor_DisplT.h"
#include "MFGP_MFA_Data_Processor_PlastT.h"

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

		/* constructors */
		MFGP_Bal_EqT 	( void );
		MFGP_Bal_EqT 	( int&, D3MeshFreeShapeFunctionT*, D3MeshFreeShapeFunctionT*, 
						GRAD_MRSSKStV*,  
						int &fTime_Step, double fdelta_t = 0.0);
		
		/* destructor */				
		~MFGP_Bal_EqT 	( void );

		void 	Initialize 	( int&, D3MeshFreeShapeFunctionT*, D3MeshFreeShapeFunctionT*, 
							GRAD_MRSSKStV*, 
							int &fTime_Step, double fdelta_t = 0.0); 

  		void 	Form_LHS_Klambda_Ku( dMatrixT &Klambda, dMatrixT &Ku); // add delta_t for dynamics
  		void 	Form_RHS_F_int( dArrayT  &F_int ); 
		void 	Form_B_List( void );  // Strain Displacement Matricies
 		void 	Form_C_List( GRAD_MRSSKStV *GRAD_MR_Plast_Mat);  // Constant List


	protected:
		dMatrixT B1_d, B3_d; 
		dMatrixT B4_lam, phi_lam; //dimension??
		//dArrayT phi_lam;
  		dMatrixT Cuu1, Cuu2, Culam1, Culam2;
  		dMatrixT Clamu1, Clamu2;
  		//double Clamlam1, Clamlam2;
  		dMatrixT Clamlam1, Clamlam2;
  		
	protected:

		MFGP_MFA_Data_Processor_DisplT Data_Pro_Displ;
		MFGP_MFA_Data_Processor_PlastT Data_Pro_Plast;

		double delta_t;
		int time_step;
		
		dSymMatrixT stress;
		dMatrixT moduli;

		int ip, n_rows_vector, n_rows_matrix, n_cols_matrix, n_sd, n_en_displ, n_en_plast, n_sd_x_n_sd, 
			n_sd_x_n_en_displ, n_sd_x_n_en_plast, Time_Integration_Scheme;
  
};


} // namespace Tahoe 
#endif /* _MFGP_BAL_EQT_H_ */

