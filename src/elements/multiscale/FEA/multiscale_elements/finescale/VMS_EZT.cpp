// $Id: VMS_EZT.cpp,v 1.4 2003-03-07 22:24:02 creigh Exp $
#include "FEA.h" 
#include "VMS.h" 

using namespace Tahoe;

VMS_EZT::VMS_EZT	( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
										int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fTime_Step,fdelta_t,Integration_Scheme);
}
 
//---------------------------------------------------------------------

void VMS_EZT::Initialize ( int &in_ip, int &in_sd, int &in_en, int Initial_Time_Step )
{
#pragma unused(in_ip)
#pragma unused(in_sd)
#pragma unused(Initial_Time_Step)
}

//---------------------------------------------------------------------

void VMS_EZT::Construct (	FEA_ShapeFunctionT &Shapes,VMF_MaterialT *EZ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
													int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
#pragma unused(fTime_Step)
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
 
	neg_1by10 = -0.1;
	//neg_1by10 = -0.0005;  // good for grad(ub) ~ grad(u)

	Data_Pro.Construct ( Shapes.dNdx );
  
  Form_A_S_Lists	(np1,n);
  Form_B_List			();

	//B.Print("B Matricies");
	//A.Print("A Matricies");

	Integral.Construct (Shapes.j, Shapes.W); 
}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_EZT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 	Ka = Integral.of( B[kB_1hat], 						B[kB_1hat] ); 
	Kb = Integral.of( B[kB_1hat], neg_1by10, 	B[kB_1hat] );
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
 	Data_Pro.grad_u ( B[kB_1hat], FEA::kNonSymmetric 	); 
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

	A[kG2]  = A[kgrad_ub];
	A[kG2] *= neg_1by10; 
	A[kG2] += A[kgrad_ua];  // G2 = grad(ua) - (1/10)grad(ub)  (i.e. Strong Form)

}


