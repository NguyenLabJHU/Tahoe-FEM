// $Id: VMS_EZ3T.cpp,v 1.1 2003-03-17 22:55:53 creigh Exp $
#include "FEA.h" 
#include "VMS.h" 

using namespace Tahoe;

VMS_EZ3T::VMS_EZ3T	(FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
						int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fTime_Step,fdelta_t,Integration_Scheme);
}

/* destructor */
//VMS_EZ3T::~VMS_EZ3T(void)
//{
//}

//---------------------------------------------------------------------

void VMS_EZ3T::Initialize (int &in_ip,int &in_sd,int &in_en, int Initial_Time_Step)
{
  n_ip = in_ip;  			// Note: Need to call Initialize() for each elmt set
  n_sd = in_sd;
	n_en = in_en;
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	C.Dimension 	(	kNUM_C_TERMS );
	A.Construct 	( kNUM_A_TERMS, 	n_ip, n_sd, n_sd );
	T4.Construct 	( kNUM_T4_TERMS, 	n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

	time_step = Initial_Time_Step;
}

//---------------------------------------------------------------------

void VMS_EZ3T::Construct (FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
#pragma unused(Integration_Scheme)

	delta_t = fdelta_t;

	Data_Pro.Construct ( Shapes.dNdx 	);
	Data_Pro.Insert_N  ( Shapes.N 		);
	Integral.Construct (Shapes.j, Shapes.W); 

  //-- The following don't allocate, that's done in Initialize()
	
  Form_C_List			( BCJ_Matl );
  Form_A_S_Lists	( np1,n );
  Form_T4_List		( );
  Form_B_List			( );

	//cout <<"C Constants: "<<C<<"\n\n"; 
	//B.Print("B Matricies");
	//A.Print("A Matricies");
	//S.Print("S Components");
	//T4.Print("T4 Matricies");

}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_EZ3T::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
	double alpha = C[kf] * -delta_t;

 	Ka  = Integral.of( B[kB_1hat], B[kB05_tau_3hat] ); 
 	Ka += Integral.of( B[kB_1hat], B[kB06_3hat]  		); 

	Ka += Integral.of( B[kB_1hat], alpha, B[kBa_S] );
	Kb  = Integral.of( B[kB_1hat], alpha, B[kBb_S] );

}

//---------------------------------------------------------------------
// F internal (F_int) dimensions here won't actually be in terms of Force

void VMS_EZ3T::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
	FEA_dVectorT G2_vec		( n_ip, n_sd_x_n_sd ); 
	Data_Pro.Reduce_Order	(	A[kG2], G2_vec 		); 

	F_int  = Integral.of	( B[kB_1hat], G2_vec	);  	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMS_EZ3T::Form_B_List (void)
{

		B.Construct (kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en);  // B = B(9,24)	

	 	Data_Pro.grad_u     	( B[kB_1hat], FEA::kNonSymmetric 	); 
	 	Data_Pro.A_grad_u_T_B	( A[kFbT], 		A[kDa_m], B[kB05_tau_3hat] 		);
	 	Data_Pro.A_grad_u_B 	( A[kDa_m], 	A[kFb],   B[kB06_3hat] 				);

		// NOTE: the term "del" is omitted before each variable for notational ease
		
	 	Data_Pro.A_grad_u_B 	( A[kA1], 		A[kFb],   B[kBa_Cb_3hat] 			); 		B[kBa_Cb_3hat] 			*= -1.0;
	 	Data_Pro.A_grad_u_B 	( A[kFbT], 		A[kFb],   B[kBb_Cb_3hat] 			);

		//-- Calculation of del ( S )
		B[kBa_S].MultAB( T4[kCC],  B[kBa_Cb_3hat] );
		B[kBb_S].MultAB( T4[kCC],  B[kBb_Cb_3hat] );

}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMS_EZ3T::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n)
{
	//---- Developer cheat: put npt in function door in-lieu-of np1 for speed

	//-----
	
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);
	A[kFai] 			= npt.Get (	VMS::kFai				);
  A[kFb]  			= npt.Get (	VMS::kFb				);	 	
	A[kFbT].Transpose  			( A[kFb] 					);
	A[kCb].MultATB     			( A[kFb], A[kFb] 	);	

	//--- LHS: del(Da) terms
	
	A[kFa_n]  = n.Get				(	VMS::kFa );
	A[kCa_n].MultATB				( A[kFa_n], 	A[kFa_n] );
	A[kDa_m].MultATBC 			(	A[kFai], 		A[kCa_n], 		A[kFai]	);
	A[kDa_m] *= 0.5;

	//--- LHS: del(sym(Eb)) terms

	A[kF_sharp].MultAB  		( A[kgrad_ub], 	A[kFb] 				);
  A[kA1].MultAB						(	A[kFbT], 			A[kF_sharp] 	);

	//=============================== RHS: Fint terms
	
	A[kDa_mp1] = 0.0;
	A[kDa_mp1].PlusIdentity(0.5); 
	
	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity (-1.0); 
	A[kEb] *= 0.5; 

	//--- Get S	
	Data_Pro.C_IJKL_E_KL 	( C[kMu], C[kLamda], A[kEb], A[kS] );	// Doesn't multiply excessive zeros 

	//-- RHS Start ******
	
  A[kG2]  = A[kS];
  A[kG2] *= (-delta_t * C[kf]); 
  A[kG2] -= A[kDa_m];
  A[kG2] += A[kDa_mp1];

	//-- RHS Finished  ******

}

//##################################################################################
//################ 4th Order Tensor Terms (i.e. T4[i] ) ########################
//##################################################################################

void VMS_EZ3T::Form_T4_List (void)  // These matricies are all 9x9 (not 3x3 like A)
{
  Data_Pro.C_IJKL (	C[kLamda],  C[kMu],					T4[kCC]			); 
}

//##################################################################################
//################ Material Constant Terms C[i] #####################################
//##################################################################################

void VMS_EZ3T::Form_C_List (VMF_MaterialT *BCJ_Matl)
{

	//-- Elasticity 	
	C[kLamda]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda 	);
	C[kMu]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu 		);

	//-- BCJ
	C[kf]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kf	 		);
}


//##################################################################################

void VMS_EZ3T::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "S" )
		tensor = A[kS];
	else if ( Name == "Eb" )
		tensor = A[kEb];
	else
		cout << " ...ERROR: VMS_EZ3T::Get() >> Unknown variable requested. \n";
}


