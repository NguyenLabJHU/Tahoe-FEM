//DEVELOPMENT

#include "FEA.h" 
#include "VMS.h" 


using namespace Tahoe;

VMS_EZT::VMS_EZT	(FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
						double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fdelta_t,Integration_Scheme);
}

//---------------------------------------------------------------------

void VMS_EZT::Construct (FEA_ShapeFunctionT &Shapes,VMF_MaterialT *EZ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
							double fdelta_t, int Integration_Scheme) 
{
#pragma unused(fdelta_t)
#pragma unused(Integration_Scheme)
#pragma unused(EZ_Matl)

	n_ip 		= np1.fVars[VMS::kGRAD_ua].IPs(); 
	n_rows	= np1.fVars[VMS::kGRAD_ua].Rows(); 
	n_cols	= np1.fVars[VMS::kGRAD_ua].Cols();
	n_en    = Shapes.dNdx.Cols();
  n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;
 
	neg_1by10 = -1.0/10.0;

	Data_Pro.Construct ( Shapes.dNdx );
  
  Form_A_S_Lists	(np1,n);
  Form_T4_List		();
  Form_B_List			();

	//B.Print("B Matricies");
	//A.Print("A Matricies");
	//T4.Print("T4 Matricies");

	Integral.Construct (Shapes.j, Shapes.W); 
}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_EZT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 /* del(grad_w)	*/ 	Ka -= Integral.of( B[kB_1hat], B[kB_IA_tau_2bar] );  	
										Kb -= Integral.of( B[kB_1hat], B[kB_IA_tau_2bar] ); 

 /* del(G2) 	  */	Ka += Integral.of( B[kB_1hat], neg_1by10, B[kB_IB1_2bar] ); 						
										Kb += Integral.of( B[kB_1hat], neg_1by10, B[kB_IB2_2bar] ); 						
										Ka += Integral.of( B[kB_1hat], neg_1by10, B[kB_IB3_2bar] ); 						
										Kb += Integral.of( B[kB_1hat], neg_1by10, B[kB_IB4_2bar] ); 						

 /* del(j)      */  Ka += Integral.of( B[kB_1hat], T4[kUI], B[kB_1hat] );
 										Kb += Integral.of( B[kB_1hat], T4[kUI], B[kB_1hat] );
}

//---------------------------------------------------------------------
// F internal (F_int) dimensions here won't actually be in terms of Force

void VMS_EZT::Form_RHS_F_int ( dArrayT &F_int) // Untested
{
	FEA_dVectorT G2_vec		( n_ip, n_sd_x_n_sd ); 
	Data_Pro.Reduce_Order	(	A[kG2], G2_vec 		); 

	F_int  = Integral.of	( B[kB_1hat], G2_vec	);  	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMS_EZT::Form_B_List (void)
{

		B.Construct (kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en);  // B = B(9,24)	

    // Generic B 
		
	 	Data_Pro.grad_u     ( B[kB_1hat], FEA::kNonSymmetric 	); 

    // del ( grad_w )
	
	 	Data_Pro.A_grad_u_T ( A[kG2], B[kB_IA_tau_2bar] ); 

		// del ( G2 )
		
	 	Data_Pro.A_grad_u 	( A[k1_minus_grad_ua], 	B[kB_IB1_2bar] ); 
	 	Data_Pro.A_grad_u 	( A[kgrad_ua], 					B[kB_IB2_2bar] ); 
	 	Data_Pro.A_grad_u 	( A[kgrad_ub], 					B[kB_IB3_2bar] ); 
	 	Data_Pro.A_grad_u 	( A[k1_minus_grad_ub], 	B[kB_IB4_2bar] ); 
		B[kB_IB2_2bar] *= -1.0;
		B[kB_IB3_2bar] *= -1.0;
		
}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

void VMS_EZT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n)
{
#pragma unused(n)

	A.Construct ( kNUM_A_TERMS, n_ip, n_sd, n_sd );	

	//-----
	
	A[kgrad_ua] 	= npt.Get (	VMS::kgrad_ua		);
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);

	A[k1_minus_grad_ua].Identity(); 
	A[k1_minus_grad_ua]  -= A[kgrad_ua];
	A[k1_minus_grad_ub].Identity(); 
	A[k1_minus_grad_ub]  -= A[kgrad_ub];

	A[kG2]  = A[kgrad_ub];
	A[kG2] *= neg_1by10; 
	A[kG2] += A[kgrad_ua];

}

//##################################################################################
//################ 4th Order Tensor Terms (i.e. T4[i] ) ########################
//##################################################################################

void VMS_EZT::Form_T4_List (void)  // These matricies are all 9x9 (not 3x3 like A)
{
	
	T4.Construct ( kNUM_T4_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

  Data_Pro.A_o_1 	( A[kG2], T4[kUI] );

}


