// $Id: VMS_BCJT.cpp,v 1.11 2003-02-03 04:40:27 paklein Exp $
#include "FEA.h" 
#include "VMS.h" 

using namespace Tahoe;

VMS_BCJT::VMS_BCJT	(FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
						double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fdelta_t,Integration_Scheme);
}

/* destructor */
//VMS_BCJT::~VMS_BCJT(void)
//{
//}

//---------------------------------------------------------------------

void VMS_BCJT::Construct (FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
							double fdelta_t, int Integration_Scheme) 
{
#pragma unused(Integration_Scheme)

	delta_t = fdelta_t;
	n_ip 		= np1.fVars[VMS::kGRAD_ua].IPs(); 
	n_rows	= np1.fVars[VMS::kGRAD_ua].Rows(); 
	n_cols	= np1.fVars[VMS::kGRAD_ua].Cols();
	n_en    = Shapes.dNdx.Cols();
  n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;


	Data_Pro.Construct ( Shapes.dNdx );
 
  Form_C_List			( BCJ_Matl );
  Form_A_S_Lists	( np1,n );
  Form_T4_List		( );
  Form_B_List			( );

	//cout <<"C Constants: "<<C<<"\n\n"; 
	//B.Print("B Matricies");
	//A.Print("A Matricies");
	//S.Print("S Components");
	//T4.Print("T4 Matricies");

	Integral.Construct (Shapes.j, Shapes.W); 
}

//---------------------------------------------------------------------
/**  Form stiffness matricies k^Alpha and k^Beta 
 *   The ON swithes in this function are for Alpha and Beta 
 *   contributions respectively */

void VMS_BCJT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 	Ka  = Integral.of( B[kB_1hat], B[kB05_tau_3hat] ); 
 	Ka += Integral.of( B[kB_1hat], B[kB06_3hat]  		); 

	Ka += Integral.of( B[kB_1hat], C[kNeg_dt_Root3by2_f], T4[kMM], B[kBa_DEV_S] );
	Kb  = Integral.of( B[kB_1hat], C[kNeg_dt_Root3by2_f], T4[kMM], B[kBb_DEV_S] );

}

//---------------------------------------------------------------------
// F internal (F_int) dimensions here won't actually be in terms of Force

void VMS_BCJT::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
	FEA_dVectorT G2_vec		( n_ip, n_sd_x_n_sd ); 
	Data_Pro.Reduce_Order	(	A[kG2], G2_vec 		); 

	F_int  = Integral.of	( B[kB_1hat], G2_vec	);  	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMS_BCJT::Form_B_List (void)
{

		B.Construct (kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en);  // B = B(9,24)	

	 	Data_Pro.grad_u     	( B[kB_1hat], FEA::kNonSymmetric 	); 
	 	Data_Pro.A_grad_u_T_B	( A[kFbT], 		A[kDa_m], B[kB05_tau_3hat] 		);
	 	Data_Pro.A_grad_u_B 	( A[kDa_m], 	A[kFb],   B[kB06_3hat] 				);

    //===================  del ( DEV( S ) ) 
	
		//--- Preparatory data
		
	 	Data_Pro.A_grad_u_B 	( A[kA1], 		A[kFb],   B[kBa_Cb_3hat] 			); 		B[kBa_Cb_3hat] 			*= -1.0;
	 	Data_Pro.A_grad_u_B 	( A[kFbT], 		A[kFb],   B[kBb_Cb_3hat] 			);
	 	Data_Pro.A_grad_u_T_B ( A[kFbT],  	A[kA1T],  B[kBa_Cb_tau_3hat] 	);		B[kBa_Cb_tau_3hat]	*= -1.0; 	
	 	Data_Pro.A_grad_u_T_B ( A[kFbT],  	A[kFb],  	B[kBb_Cb_tau_3hat] 	);
	  Data_Pro.A_grad_u_B 	( A[kA2],   	A[kA3],   B[kBa_Cbi_3hat] 		);
	 	Data_Pro.A_grad_u_B 	( A[kA3T],  	A[kA3],	  B[kBb_Cbi_3hat] 		);		B[kBb_Cbi_3hat]			*= -1.0; 	
	 	Data_Pro.A_grad_u_T_B ( A[kA3T],  	A[kA2T],  B[kBa_Cbi_tau_3hat] );
	 	Data_Pro.A_grad_u_T_B ( A[kA3T],  	A[kA3],  	B[kBb_Cbi_tau_3hat] );		B[kBb_Cbi_tau_3hat]	*= -1.0; 	

		B[kBa_S].MultAB( T4[kCC],  B[kBa_Cb_3hat] );
		B[kBb_S].MultAB( T4[kCC],  B[kBb_Cb_3hat] );

		//--- Calculation of Ba_DEV_S
		
		T4[kT4_Temp0]  = T4[kCbi_o_Cb];		
		T4[kT4_Temp0] *= -C[k1by3]; 
		T4[kT4_Temp0].PlusIdentity(); 
		B[kBa_DEV_S].MultAB( 	T4[kT4_Temp0], B[kBa_S] );

	  B[kB_Temp0].SumOf  ( 	B[kBa_Cb_3hat], B[kBa_Cb_tau_3hat]	) ;		
		B[kB_Temp1].MultAB ( 	T4[kCbi_o_S], B[kB_Temp0]  );
		B[kB_Temp1]  	*= 		  C[k1by3]; 
		B[kBa_DEV_S] 	-=			B[kB_Temp1]; 

	  B[kB_Temp0].SumOf( 		B[kBa_Cbi_3hat], B[kBa_Cbi_tau_3hat]	) ;	
		B[kB_Temp0] 	*=  		S[kCb_i_S]; 
		B[kB_Temp0] 	*= 		 	C[k1by3]; 
		B[kBa_DEV_S] 	-=			B[kB_Temp0]; 

		//--- Calculation of Bb_DEV_S
		
		B[kBb_DEV_S].MultAB	( T4[kT4_Temp0], B[kBb_S] );

	  B[kB_Temp0].SumOf  ( 	B[kBb_Cb_3hat], B[kBb_Cb_tau_3hat]	) ;		
		B[kB_Temp1].MultAB ( 	T4[kCbi_o_S], B[kB_Temp0]  );
		B[kB_Temp1]  	*= 		  C[k1by3]; 
		B[kBb_DEV_S] 	-=			B[kB_Temp1]; 

	  B[kB_Temp0].SumOf( 		B[kBb_Cbi_3hat], B[kBb_Cbi_tau_3hat]	) ;	
		B[kB_Temp0] 	*=  		S[kCb_i_S]; 
		B[kB_Temp0] 	*= 		 	C[k1by3]; 
		B[kBb_DEV_S] 	-=			B[kB_Temp0]; 

}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMS_BCJT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n)
{
	//---- Developer cheat: put npt in function door in-lieu-of np1 for speed

	A.Construct ( kNUM_A_TERMS, n_ip, n_sd, n_sd );
	S.Construct ( kNUM_S_TERMS, n_ip );

	//-----
	
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);
	A[kFai] 			= npt.Get (	VMS::kFai				);
  A[kFb]  			= npt.Get (	VMS::kFb				);	 	
	A[kFbT].Transpose  			( A[kFb] 					);

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
	
	A[kCb].MultATB     			( A[kFb], A[kFb] );	
	A[kCbi].Inverse    			( A[kCb] ); 

	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity (-1.0); 
	A[kEb] *= 0.5; 

	//--- Get S	
	Data_Pro.C_IJKL_E_KL 	( C[kMu], C[kLamda], A[kEb], A[kS] );	// Doesn't multiply excessive zeros 

	/*A[kgrad_ub].Print("grad_ub");
	A[kEb].Print("Eb");
	A[kS].Print("S"); */

	//--- Get DEV(S) : Put this as a method in DataPro ::DEV(Fb,S) ( use MultABCT(Fb,dev(S),Fb) )
		S[kCb_i_S].Double_Dot ( A[kCb], A[kS] ); 
		A[kDEV_S]  = 	A[kCbi];
		A[kDEV_S] *= 	S[kCb_i_S]; 
		A[kDEV_S] *= -C[k1by3]; 
		A[kDEV_S] += 	A[kS]; 

	//--- Get ||DEV(S)||, N
  A[kDEV_S].Mag_and_Dir ( S[kMag_DEV_S], A[kN] );

	//--- Get Beta
	S[kBeta]  = S[kMag_DEV_S]; 
	S[kBeta] *= C[kRoot3by2byV]; 
	//S[kBeta].Print(">>### kBeta ###<<");

	//--- Get Sinh(Beta) 
	S[kSinh_Beta].Sinh( S[kBeta] ); 
	S[kCosh_Beta].Cosh( S[kBeta] ); 
	//S[kSinh_Beta].Print(">>### kSinh_Beta ###<<");
	//S[kCosh_Beta].Print(">>### kCosh_Beta ###<<");

	int sinh_ON=1;
	//--- Get G2
  A[kG2]  = A[kN];
	if (sinh_ON)
  	A[kG2] *= S[kSinh_Beta]; 
	else
  	A[kG2] *= S[kBeta]; 
  A[kG2] *= C[kNeg_dt_Root3by2_f]; 
  A[kG2] -= A[kDa_m];
  A[kG2] += A[kDa_mp1];

	//=============================== Misc LHS
	
	//--- Get Cosh(Beta), Beta2
 	S[kBeta2]  = C[kRoot3by2byV]; 
	if (sinh_ON)
		S[kBeta2]  *= S[kCosh_Beta];

}

//##################################################################################
//################ 4th Order Tensor Terms (i.e. T4[i] ) ########################
//##################################################################################

void VMS_BCJT::Form_T4_List (void)  // These matricies are all 9x9 (not 3x3 like A)
{
	
	T4.Construct ( kNUM_T4_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

  Data_Pro.C_IJKL (	C[kLamda],  C[kMu],					T4[kCC]			); 

	A[kA_Temp0] = A[kN];
	Data_Pro.A_o_B 					( A[kCbi], 	A[kS], 				T4[kCbi_o_S] 	);
	Data_Pro.A_o_B 					( A[kCbi], 	A[kCb], 			T4[kCbi_o_Cb] );
	Data_Pro.A_o_B 					( A[kN], 		A[kA_Temp0], 	T4[kN_o_N] 		);
 	Data_Pro.II_minus_A_o_B ( A[kN],		A[kA_Temp0],	T4[kPP] 			); // 4th Order Projector 

	//-- Build MM
  T4[kMM]  = T4[kPP];
  T4[kMM] *= S[kMag_DEV_S];
	T4[kT4_Temp0]  = T4[kN_o_N];
	T4[kT4_Temp0] *= S[kBeta2]; 
  T4[kMM] += T4[kT4_Temp0];

}

//##################################################################################
//################ Material Constant Terms C[i] #####################################
//##################################################################################

void VMS_BCJT::Form_C_List (VMF_MaterialT *BCJ_Matl)
{

	C.Dimension (kNUM_C_TERMS);
	
	C[kLamda]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda 	);
	C[kMu]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu 		);
	C[kf]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kf	 		);
	C[kV]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kV	 		);
	C[kY]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kY	 		);

	C[kRoot3by2]						= sqrt(1.5);
	C[kNeg_dt_Root3by2_f]  	= -1.0 * delta_t * C[kRoot3by2] * C[kf]; 
	C[kRoot3by2byV]					= C[kRoot3by2] / C[kV];
	C[k1byV]								=	1.0 /  C[kV];
	C[kYbyV]								=	C[kY] /  C[kV];

	C[k1by3] = 1.0/3.0;

}


