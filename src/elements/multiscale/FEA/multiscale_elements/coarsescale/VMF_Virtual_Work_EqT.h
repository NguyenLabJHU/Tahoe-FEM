// $Id: VMF_Virtual_Work_EqT.h,v 1.5 2003-02-03 04:40:26 paklein Exp $
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
						   	kS,
						   	kSigma,
						   	kA_Temp0,
            		kNUM_A_TERMS };  

    	enum S_T { 
								kJ,
	             	kNUM_S_TERMS };  
						
    	enum T4_T { 
								kd_1bar,
						   	kT4_Temp0,
	             	kNUM_T4_TERMS };  
						

 
		//--------------------------------------------------------------
		
		VMF_Virtual_Work_EqT 	( void ) { } 

		VMF_Virtual_Work_EqT 	( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
														int Integration_Scheme=FEA::kBackward_Euler);

		void 	Construct 			( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
														int Integration_Scheme=FEA::kBackward_Euler); 

  	void 	Form_LHS_Ka_Kb	(	dMatrixT &Ka, dMatrixT &Kb, double delta_t =0.0	); // add delta_t for dynamics
  	void 	Form_RHS_F_int	(	dArrayT  &F_int, double delta_t =0.0	); 
		void 	Form_B_List 		( void );  // Strain Displacement Matricies
		void 	Form_A_S_Lists 	( VMS_VariableT &np1, VMS_VariableT &n,int Integration_Scheme=FEA::kBackward_Euler ); // BCDE ---> A 

		void 	Sigma ( FEA_dMatrixT &sigma) 		{ sigma = A[kSigma]; } 

	protected:

  	FEA_dMatrix_ArrayT B; 
  	FEA_dMatrix_ArrayT A; 
  	FEA_dMatrix_ArrayT T4; 
  	FEA_dScalar_ArrayT S; 

	protected:

		FEA_IntegrationT 		Integral;
		FEA_Data_ProcessorT Data_Pro; 

		double lamda,mu;
		int n_ip, n_rows, n_cols, n_sd, n_en, n_sd_x_n_sd, n_sd_x_n_en;
  
};

} // namespace Tahoe 
#endif /* _VMS_BCJT_H_ */

