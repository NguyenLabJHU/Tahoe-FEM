// $Id: VMS_EZT.cpp,v 1.5 2003-03-17 22:05:32 creigh Exp $
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
#pragma unused(Initial_Time_Step)
	n_ip = in_ip;
	n_sd = in_sd;
	n_en = in_en;
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	A.Construct ( kNUM_A_TERMS, n_ip, n_sd, n_sd );	
	B.Construct ( kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en);  // B = B(9,24)	
	C.Dimension ( kNUM_C_TERMS );

}

//---------------------------------------------------------------------

void VMS_EZT::Construct (	FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
													int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
#pragma unused(fTime_Step)
#pragma unused(fdelta_t)
#pragma unused(Integration_Scheme)

 
	//C[kNeg_Alpha] = -0.0005;  // good for grad(ub) ~ grad(u)

	Data_Pro.Construct ( Shapes.dNdx );
  
  Form_C_List			( BCJ_Matl 	);
  Form_A_S_Lists	( np1,n 		);
  Form_B_List			();

	Integral.Construct (Shapes.j, Shapes.W); 
}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_EZT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 	Ka = Integral.of( B[kB_1hat], B[kB_1hat] ); 
	Kb = Integral.of( B[kB_1hat], C[kNeg_Alpha], 	B[kB_1hat] );
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
 	Data_Pro.grad_u ( B[kB_1hat], FEA::kNonSymmetric 	); 
}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

void VMS_EZT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n)
{
#pragma unused(n)

	A[kgrad_ua] 	= npt.Get (	VMS::kgrad_ua		);
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);

	A[kG2]  = A[kgrad_ub];
	A[kG2] *= C[kNeg_Alpha]; 
	A[kG2] += A[kgrad_ua];  // G2 = grad(ua) - alpha.grad(ub)  (i.e. Strong Form)

}

//##################################################################################
//################ C_TERMS #########################################################
//##################################################################################

void VMS_EZT::Form_C_List (VMF_MaterialT *BCJ_Matl)
{
	C[kNeg_Alpha] = BCJ_Matl -> Retrieve ( BCJ_MatlT::kf	);
	C[kNeg_Alpha] *= -1.0; 
}

//##################################################################################

void VMS_EZT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "grad_ua" )
		tensor = A[kgrad_ua];
	else if ( Name == "grad_ub" )
		tensor = A[kgrad_ub];
	else
		cout << " ...ERROR: VMS_EZ3T::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}

