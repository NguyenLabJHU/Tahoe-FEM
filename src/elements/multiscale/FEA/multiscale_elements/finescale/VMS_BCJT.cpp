// $Id: VMS_BCJT.cpp,v 1.12 2003-03-07 22:24:02 creigh Exp $
#include "FEA.h" 
#include "VMS.h" 

using namespace Tahoe;

VMS_BCJT::VMS_BCJT	(FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
						int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,BCJ_Matl,np1,n,fTime_Step,fdelta_t,Integration_Scheme);
}

/* destructor */
//VMS_BCJT::~VMS_BCJT(void)
//{
//}

//---------------------------------------------------------------------

void VMS_BCJT::Initialize (int &in_ip,int &in_sd,int &in_en, int Initial_Time_Step)
{
  n_ip = in_ip;  			// Note: Need to call Initialize() for each elmt set
  n_sd = in_sd;
	n_en = in_en;
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	C.Dimension 	(	kNUM_C_TERMS );
	S.Construct 	( kNUM_S_TERMS, 	n_ip 	);
	A.Construct 	( kNUM_A_TERMS, 	n_ip, n_sd, n_sd );
	T4.Construct 	( kNUM_T4_TERMS, 	n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

	NP.Construct	( kNUM_NP_TERMS, n_en, n_sd, n_sd );  // Special one (not FEA at ip) rather values at Node Point (NP)
	NP[kIdentity].Identity( );

	S[kIV_Alpha_n] 	= 0.0; 		//-- Initialize for flow rule 
	S[kIV_Alpha] 		= 0.0;   	//-- Why not?

	time_step = Initial_Time_Step;
}

//---------------------------------------------------------------------

void VMS_BCJT::Construct (FEA_ShapeFunctionT &Shapes,VMF_MaterialT *BCJ_Matl, VMS_VariableT &np1, VMS_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
#pragma unused(Integration_Scheme)


	if ( fTime_Step != time_step) { 	// New time step
		S[kIV_Alpha_n] = S[kIV_Alpha];  // Eventually do this with all the kVar_n and make method Next_Step()
		time_step = fTime_Step;
	}

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

void VMS_BCJT::Form_LHS_Ka_Kb ( dMatrixT &Ka, dMatrixT &Kb )
{
 	Ka  = Integral.of( B[kB_1hat], B[kB05_tau_3hat] ); 
 	Ka += Integral.of( B[kB_1hat], B[kB06_3hat]  		); 

	Ka += Integral.of( B[kB_1hat], C[kNeg_dt_Root3by2_f], T4[kMM], B[kBa_DEV_H] );
	Kb  = Integral.of( B[kB_1hat], C[kNeg_dt_Root3by2_f], T4[kMM], B[kBb_DEV_H] );

 	if 	(Iso_Hard_Type != kNo_Iso_Hard) 
		Ka += Integral.of( B[kB_1hat], S[kBetaK], T4[kN_o_Ne], B[kBa_Kappa] );
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

    //===================  del ( DEV( H ) ) 
		
		// NOTE: the term "del" is omitted before each variable for notational ease
		
		//--- Preparatory data
		
	 	Data_Pro.A_grad_u_B 	( A[kA1], 		A[kFb],   B[kBa_Cb_3hat] 			); 		B[kBa_Cb_3hat] 			*= -1.0;
	 	Data_Pro.A_grad_u_B 	( A[kFbT], 		A[kFb],   B[kBb_Cb_3hat] 			);
	 	Data_Pro.A_grad_u_T_B ( A[kFbT],  	A[kA1T],  B[kBa_Cb_tau_3hat] 	);		B[kBa_Cb_tau_3hat]	*= -1.0; 	
	 	Data_Pro.A_grad_u_T_B ( A[kFbT],  	A[kFb],  	B[kBb_Cb_tau_3hat] 	);
	  Data_Pro.A_grad_u_B 	( A[kA2],   	A[kA3],   B[kBa_Cbi_3hat] 		);
	 	Data_Pro.A_grad_u_B 	( A[kA3T],  	A[kA3],	  B[kBb_Cbi_3hat] 		);		B[kBb_Cbi_3hat]			*= -1.0; 	
	 	Data_Pro.A_grad_u_T_B ( A[kA3T],  	A[kA2T],  B[kBa_Cbi_tau_3hat] );
	 	Data_Pro.A_grad_u_T_B ( A[kA3T],  	A[kA3],  	B[kBb_Cbi_tau_3hat] );		B[kBb_Cbi_tau_3hat]	*= -1.0; 	

		//-- Calculation of del ( S )
		B[kBa_S].MultAB( T4[kCC],  B[kBa_Cb_3hat] );
		B[kBb_S].MultAB( T4[kCC],  B[kBb_Cb_3hat] );

		//-- Calculation of del ( Zeta )
		Form_del_Zeta_B ( ); // Calculates B[kBa_Z],  B[kBb_Z]	

		B[kBa_H].DiffOf (	B[kBa_S], B[kBa_Z] );
		B[kBb_H].DiffOf (	B[kBb_S], B[kBb_Z] );

		//--- Calculation of Ba_DEV_H : del (DEV(H)) = (1-Cbi_o_Cb):del(H) - 1/3(CbiT:del(Cb)) - 1/3.BetaTC.del(Cbi)
		
		T4[kT4_Temp0]  = T4[kCbi_o_Cb];		
		T4[kT4_Temp0] *= -C[k1by3]; 
		T4[kT4_Temp0].PlusIdentity(); 
		B[kBa_DEV_H].MultAB( 	T4[kT4_Temp0], B[kBa_H] );

	  B[kB_Temp0].SumOf  ( 	B[kBa_Cb_3hat], B[kBa_Cb_tau_3hat]	) ;		
		B[kB_Temp1].MultAB ( 	T4[kCbi_o_H], B[kB_Temp0]  );
		B[kB_Temp1]  	*= 		  C[k1by3]; 
		B[kBa_DEV_H] 	-=			B[kB_Temp1]; 

	  B[kB_Temp0].SumOf( 		B[kBa_Cbi_3hat], B[kBa_Cbi_tau_3hat]	) ;	
		B[kB_Temp0] 	*=  		S[kCb_i_H]; 
		B[kB_Temp0] 	*= 		 	C[k1by3]; 
		B[kBa_DEV_H] 	-=			B[kB_Temp0]; 

		//--- Calculation of Bb_DEV_H
		
		B[kBb_DEV_H].MultAB	( T4[kT4_Temp0], B[kBb_H] );

	  B[kB_Temp0].SumOf  ( 	B[kBb_Cb_3hat], B[kBb_Cb_tau_3hat]	) ;		
		B[kB_Temp1].MultAB ( 	T4[kCbi_o_H], B[kB_Temp0]  );
		B[kB_Temp1]  	*= 		  C[k1by3]; 
		B[kBb_DEV_H] 	-=			B[kB_Temp1]; 

	  B[kB_Temp0].SumOf( 		B[kBb_Cbi_3hat], B[kBb_Cbi_tau_3hat]	) ;	
		B[kB_Temp0] 	*=  		S[kCb_i_H]; 
		B[kB_Temp0] 	*= 		 	C[k1by3]; 
		B[kBb_DEV_H] 	-=			B[kB_Temp0]; 

    //===================  del ( DEV( H ) ) [END]
	
		//--- Calculation of Ba_Kappa_bar ( Iso Hardening )

 		if 	(Iso_Hard_Type != kNo_Iso_Hard) {

	 		Data_Pro.A_grad_u_B 	( A[kFaT], 	A[kF],   B[kBa_Kappa_3hat] 			);
	 		Data_Pro.A_grad_u_T_B ( A[kFT],  	A[kFa],  B[kBa_Kappa_tau_3hat] 	);	
			B[kBa_Kappa].SumOf 		( B[kBa_Kappa_tau_3hat], B[kBa_Kappa_3hat] 	);

		}
		
}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMS_BCJT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n)
{
	//---- Developer cheat: put npt in function door in-lieu-of np1 for speed

	//-----
	
	A[kF] 				= npt.Get (	VMS::kF					);
	A[kFa] 				= npt.Get (	VMS::kFa				);
	A[kgrad_ub] 	= npt.Get (	VMS::kgrad_ub		);
	A[kFai] 			= npt.Get (	VMS::kFai				);
  A[kFb]  			= npt.Get (	VMS::kFb				);	 	
  A[kFbi]  			= npt.Get (	VMS::kFbi				);	 	
	A[kFbT].Transpose  			( A[kFb] 					);
	A[kFaT].Transpose  			( A[kFa] 					);
	A[kCa].MultATB     			( A[kFa], A[kFa] 	);	
	A[kCb].MultATB     			( A[kFb], A[kFb] 	);	
	A[kCbi].Inverse    			( A[kCb] ); 
	S[kJb].Determinant 			( A[kFb] );


	//--- LHS: del(Da) terms
	
	A[kFa_n]  = n.Get				(	VMS::kFa );
	A[kCa_n].MultATB				( A[kFa_n], 	A[kFa_n] );
	A[kDa_m].MultATBC 			(	A[kFai], 		A[kCa_n], 		A[kFai]	);
	A[kDa_m] *= 0.5;

	//--- LHS: del(sym(Eb)) terms

	A[kF_sharp].MultAB  		( A[kgrad_ub], 	A[kFb] 				);
  A[kA1].MultAB						(	A[kFbT], 			A[kF_sharp] 	);
  A[kA2].MultAB						(	A[kCbi], 			A[kA1] 				);
  A[kA3].MultAB						(	A[kFb],				A[kCbi] 			);
	A[kA1T].Transpose 			( A[kA1] ); 
	A[kA2T].Transpose 			( A[kA2] ); 
	A[kA3T].Transpose 			( A[kA3] ); 

	//--- LHS: del( DEV(Zeta) ) terms
	
	A[kF_sharp_T].Transpose ( A[kF_sharp] );

	//=============================== RHS: Fint terms
	
	A[kDa_mp1] = 0.0;
	A[kDa_mp1].PlusIdentity(0.5); 
	
	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity (-1.0); 
	A[kEb] *= 0.5; 

	//--- Get Zeta	

 	if 	( Back_Stress_Type != kNo_Back_Stress ) //-- Apply Back Stress 
		Get_Back_Stress ( ); //-- Builds  A[kZeta]
	else
		A[kZeta] = 0.0;
	
	//--- Get S	
	Data_Pro.C_IJKL_E_KL 	( C[kMu], C[kLamda], A[kEb], A[kS] );	// Doesn't multiply excessive zeros 

	/* A[kgrad_ub].Print("grad_ub");
	A[kEb].Print("Eb");
	A[kS].Print("S"); */

	A[kH_bar].SumOf ( A[kS], A[kZeta] );

	//--- Get DEV(S) : Put this as a method in DataPro ::DEV(Fb,S) ( use MultABCT(Fb,dev(S),Fb) )
	S[kCb_i_H].Double_Dot ( A[kCb], A[kH_bar] ); 
	A[kDEV_H]  = 	A[kCbi];
	A[kDEV_H] *= 	S[kCb_i_H]; 
	A[kDEV_H] *= -C[k1by3]; 
	A[kDEV_H] += 	A[kH_bar]; 

	//--- Get ||DEV(H)||, N
  A[kDEV_H].Mag_and_Dir ( S[kMag_DEV_H], A[kN] );

	//--- Get Beta
	S[kBeta]  = S[kMag_DEV_H];
	S[kBeta] *= C[kRoot3by2]; 
 	if 	( Iso_Hard_Type != kNo_Iso_Hard ) { //-- Apply Iso Hard to RHS
		Get_Iso_Hard_Kappa_bar ( ); 
		S[kBeta] -= S[kKappa_bar]; 
		S[kBeta] -= C[kY]; 
	}
	S[kBeta] *= C[k1byV];

	//S[kBeta].Print(">>### kBeta ###<<");

	//--- Get Sinh(Beta) 
	S[kMacaulay_Sinh_Beta].Sinh( S[kBeta] ); 
	S[kMacaulay_Sinh_Beta].Macaulay ( );

	S[kBeta2].Cosh( S[kBeta] ); 
 	S[kBeta2] *= C[kRoot3by2byV]; 

	//--- Include Hardening Coef. Kappa_Bar
 	if 	( Iso_Hard_Type != kNo_Iso_Hard ) {
		S[kBetaK] 		= S[kBeta2]; 
		S[kBetaK] 	 *= ( delta_t * C[kf] * C[k1by3] * C[kKappa] ); //-- For LHS
	}	

	//S[kMacaulay_Sinh_Beta].Print(">>### kMacaulay_Sinh_Beta ###<<");
	//S[kBeta2].Print(">>### kcosh(beta)*root(3/2)/V ###<<");

	//-- RHS Start ******
	
  A[kG2]  = A[kN];
  A[kG2] *= S[kMacaulay_Sinh_Beta]; 
  A[kG2] *= C[kNeg_dt_Root3by2_f]; 
  A[kG2] -= A[kDa_m];
  A[kG2] += A[kDa_mp1];

	//-- RHS Finished  ******
	
}

//##################################################################################
//################ 4th Order Tensor Terms (i.e. T4[i] ) ########################
//##################################################################################

void VMS_BCJT::Form_T4_List (void)  // These matricies are all 9x9 (not 3x3 like A)
{
	
  Data_Pro.C_IJKL (	C[kLamda],  C[kMu],					T4[kCC]			); 

	A[kA_Temp0] = A[kN];
	Data_Pro.A_o_B 					( A[kCbi], 			A[kH_bar], 		T4[kCbi_o_H] 	);
	Data_Pro.A_o_B 					( A[kCbi], 			A[kCb], 			T4[kCbi_o_Cb] );
	Data_Pro.A_o_B 					( A[kA_Temp0], 	A[kN], 				T4[kN_o_N] 		);  
 	Data_Pro.II_minus_A_o_B ( A[kA_Temp0],	A[kN],				T4[kPP] 			); // 4th Order Projector 

 	if 	(Iso_Hard_Type != kNo_Iso_Hard) 
		Data_Pro.A_o_B 				( A[kA_Temp0], 	A[kNe], 	T4[kN_o_Ne] 	);  

 	if 	( Back_Stress_Type != kNo_Back_Stress ) { 
		Data_Pro.A_o_B 				( A[kZeta], 	A[kF_sharp_T], 	T4[kZ_o_F_sharp_T] 	);  
		Data_Pro.A_o_1 				( A[kZeta], 	T4[kZ_o_1] 	);  
	}

	//-- Build MM
	T4[kT4_Temp0]  = T4[kN_o_N];
	T4[kT4_Temp0] *= S[kBeta2]; 
  T4[kMM]  = T4[kPP];
  T4[kMM] *= S[kMacaulay_Sinh_Beta];
  T4[kMM] += T4[kT4_Temp0];

}
//kMag_DEV_H_PP
//##################################################################################
//################ Material Constant Terms C[i] #####################################
//##################################################################################

void VMS_BCJT::Form_C_List (VMF_MaterialT *BCJ_Matl)
{

	//-- Elasticity 	
	C[kLamda]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kLamda 	);
	C[kMu]    	 	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kMu 		);

	//-- BCJ
	C[kf]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kf	 		);
	C[kV]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kV	 		);
	C[kY]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kY	 		);

	//-- Iso Hard
	C[kKappa]    	= BCJ_Matl -> Retrieve ( BCJ_MatlT::kPlastic_Modulus_K );
	C[kH]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kH	 		);

	//-- Kine Hard (Back Stress)
	C[kc]    			= BCJ_Matl -> Retrieve ( BCJ_MatlT::kc_zeta );
	C[kl]    	 		= BCJ_Matl -> Retrieve ( BCJ_MatlT::kl );

	//-- Combinations
	C[kRoot3by2]						= sqrt(1.5);
	C[kNeg_dt_Root3by2_f]  	= -1.0 * delta_t * C[kRoot3by2] * C[kf]; 
	C[kRoot3by2byV]					= C[kRoot3by2] / C[kV];
	C[k1byV]								=	1.0 /  C[kV];
	C[kYbyV]								=	C[kY] /  C[kV];
	C[k1by3] 								= 1.0/3.0;
	C[k2by3] 								= 2.0/3.0;
	C[kRoot2by3]						= sqrt( C[k2by3] ); 
	C[kMu_c_l]							= C[kMu] * C[kc] * C[kl];

}


//##################################################################################
//################ Isotropic Hardening #############################################
//##################################################################################

void VMS_BCJT::Get_Iso_Hard_Kappa_bar ( void )
{
		A[kEa]  = A[kCa]; 
		A[kEa].PlusIdentity (-1.0); 
		A[kEa] *= 0.5; 

		A[kEa_n]  = A[kCa_n]; 
		A[kEa_n].PlusIdentity (-1.0); 
		A[kEa_n] *= 0.5; 

		A[kEa_dot].DiffOf ( A[kEa], A[kEa_n] );
		A[kEa_dot] /= delta_t; 
  	A[kEa_dot].Mag_and_Dir 		( S[kMag_Ea_dot], 	A[kNe] );

		if ( Iso_Hard_Type == kMethod1 || Iso_Hard_Type == kMethod2 ) {

			//-- Derive Kappa Bar Term 1
			S[kKappa_bar]		 = S[kIV_Alpha_n];    
			S[kKappa_bar]		*= ( C[kKappa]*C[kRoot2by3] );    

			//-- Derive Kappa Bar Term 2 
  		S[kS_Temp0]  = S[kMag_Ea_dot]; 
  		S[kS_Temp0] *= ( C[kRoot2by3] * C[kKappa] * C[kH] * delta_t ); 

			//-- Derive Kappa Bar: Sum Term 1 and Term 2 
  		S[kKappa_bar]  += S[kS_Temp0];  

		}
		else
			cout << " ...ERROR >> VMS_BCJT::Get_Iso_Hard_Kappa() : Bad Iso_Hard_Type \n";

}

//##################################################################################
//################ Kinematic Hardening ( Back Stress ) #############################
//##################################################################################

//##################################################### Back Stress (RHS) 

void VMS_BCJT::Get_Back_Stress ( void )
{
	if ( Back_Stress_Type == kSteinmann ) {

		int n_ed = n_sd_x_n_sd;
		int n_ed_x_n_en = n_ed*n_en; 
																																				// -- 2D --
		FEA_dMatrixT 	me 					( n_ip, n_ed_x_n_en, n_ed_x_n_en );  			// (16x16)
		FEA_dMatrixT 	B_tilde 		( n_ip, n_ed, n_ed_x_n_en ); 							// (4x16)
		FEA_dVectorT 	grad_ub_vec ( n_ip, n_sd_x_n_sd ); 										// (4x1)

		dMatrixT 			M						( n_ed_x_n_en, n_ed_x_n_en ); 						// (16x16)
		dArrayT				fe					( n_ed_x_n_en ); 													// (16x1)
		dArrayT				sE_vec			( n_sd_x_n_sd * n_en ); 									// (16x1)

		/*sE_mat.Dimension ( n_en ); // sE_mat is a member variable
		for (int a=0; a<n_en; a++)
			sE_mat[a].Dimension ( n_sd ); */

		Data_Pro.Reduce_Order	(	A[kgrad_ub], grad_ub_vec ); 
		Data_Pro.Mass 				( n_ed, me );
		Data_Pro.Mass_B 			( n_ed, B_tilde );

 		M  = Integral.of  ( me );
		fe = Integral.of	( B_tilde, grad_ub_vec	);  	// (16x4)(4x1) = (16x1)

		M.Inverse ( );

		M.Multx ( fe,sE_vec ); // (16x16)(16x1) = (16x1)

		Data_Pro.Element_Nodal_Values_Expand_Order ( sE_vec , NP[ksE_mat] );	// (16x1) --> (4x2x2)

		//-- Static Condensation follows:
	
		Data_Pro.Curl ( NP[ksE_mat], A[kcurl_sE] );  // Recall Curl ( 1 - sE ) = Curl ( 1 ) - Curl (sE) = -Curl(sE);

		// Note: two negatives cancel each other: Zeta = - c*mu*l*Jb*Fbi*(Curl( -sE )^T) 

		A[kZeta].MultABT ( A[kFbi], A[kcurl_sE] );
		A[kZeta] *= S[kJb]; 
		A[kZeta] *= C[kMu_c_l]; 

		//A[kZeta].Print( "Zeta" );

	}
	else
		cout << " VMS_BCJT::Get_Back_Stress() >> Unknown Back_Stress_Type \n";
	
}

//##################################################### Back Stress (LHS) [ Linearized ]
// Note: all C[],S[],A[],and T4[] terms already formed

void VMS_BCJT::Form_del_Zeta_B ( void ) // Calculates B[kBa_Z],  B[kBb_Z]	
{

 	if 	( Back_Stress_Type == kNo_Back_Stress ) {
		B[kBa_Z] = 0.0;
		B[kBb_Z] = 0.0;
		return;
	}	

	bool bDEL_sE=1, bDEL_grad_ub=1;

	//-- Mu*c*l*del(Jb)*Fb^-1*curl(sE)^T
			
	B[kBa_Zeta_Jb].MultAB ( T4[kZ_o_F_sharp_T], 	B[kB_1hat] );  	// Needs neg sign
	B[kBb_Zeta_Jb].MultAB ( T4[kZ_o_1], 					B[kB_1hat] );		// Note: Zeta = Mu_c_l_Fbi_(curl(sE))^T

	//-- Mu*c*l*Jb*del(Fb^-1)*curl(sE)^T

	A[kA4].MultAB ( A[kFbi], A[kF_sharp] );	

	A[kMu_c_l_Jb_curl_sE_T].Transpose(  A[kcurl_sE] );
	A[kMu_c_l_Jb_curl_sE_T] *= S[kJb]; 
	A[kMu_c_l_Jb_curl_sE_T] *= C[kMu_c_l]; 

	Data_Pro.A_grad_u_B ( A[kA4], 	A[kMu_c_l_Jb_curl_sE_T],   B[kBa_Zeta_Fbi] );
	Data_Pro.A_grad_u_B ( A[kFbi], 	A[kMu_c_l_Jb_curl_sE_T],   B[kBb_Zeta_Fbi] ); // Needs neg sign

	if ( bDEL_sE == bDEL_grad_ub ) { // sE behaves as grad_ub, likewise del(sE) will behave as del(grad_ub)

		//-- Mu*c*l*Jb*Fb^-1*del(curl(sE)^T)
	
		NP[kFbi_hat_mat]  = NP[kIdentity];	// Fbi_hat = 1 - sE
		NP[kFbi_hat_mat] -= NP[ksE_mat];

		int i,j,s,r,m,p;
		double E_rmp;

		FEA_dScalarT Fip_sEjsm_Ermp ( n_ip );  // F,sE,E -->  Fbi, script Epsilon, permutation sybol
		FEA_dScalarT sE_jsm 				( n_ip );  // d(sE(j,s))/dm 

		FEA_dScalarT Fip_hFjsm_Ermp ( n_ip );  // F,hF,E -->  Fbi, Fbi_hat, permutation sybol
		FEA_dScalarT Fbi_hat_jsm		( n_ip );  // d(Fbi_hat(j,s))/dm 

		FEA_dScalarT Fip_sEjms_Emrp ( n_ip );  // F,sE,E -->  Fbi, script Epsilon, permutation sybol
		FEA_dScalarT sE_jms 				( n_ip );  // d(sE(j,m))/ds 

		for (i=0; i<n_sd; i++)		// Form 4th Order Tensors
			for (j=0; j<n_sd; j++)
				for (s=0; s<n_sd; s++)
					for (r=0; r<n_sd; r++) {

						T4[kAAe_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) = 0.0; 
						T4[kAAf_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) = 0.0; 
						T4[kEE_2bar] ( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) = 0.0; 

						//-- Loop over dummy indicies m,p --> summation forms AA_2bar 4th order tensors
					
						for (m=0; m<n_sd; m++) {
							for (p=0; p<n_sd; p++) {

								E_rmp = Data_Pro.e_ijk[r][m][p];  // Permutation symbol

								Data_Pro.Grad_ij ( NP[ksE_mat],     	j,s,m, sE_jsm 			); 
								Data_Pro.Grad_ij ( NP[kFbi_hat_mat],	j,s,m, Fbi_hat_jsm 	);
								Data_Pro.Grad_ij ( NP[ksE_mat],     	j,m,s, sE_jms 			); 

								//-- AAe_2bar 
								Fip_sEjsm_Ermp  = A[kFbi](i,p); 
								Fip_sEjsm_Ermp *= sE_jsm;
								Fip_sEjsm_Ermp *= E_rmp; 
								T4[kAAe_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) +=  Fip_sEjsm_Ermp;

								//-- AAf_2bar 
								Fip_hFjsm_Ermp  = A[kFbi](i,p); 
								Fip_hFjsm_Ermp *= Fbi_hat_jsm;
								Fip_hFjsm_Ermp *=  E_rmp; 
								T4[kAAf_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) +=  Fip_hFjsm_Ermp;

								//-- EE_2bar 
								Fip_sEjms_Emrp  = A[kFbi](i,p); 
								Fip_sEjms_Emrp *= sE_jms;
								Fip_sEjms_Emrp *=  -E_rmp; // Recall: E_mrp = -E_rmp
								T4[kEE_2bar]( Data_Pro.Map(i,j) , Data_Pro.Map(s,r) ) +=  Fip_sEjms_Emrp;

					  	}	
						}
					}
	
		S[kS_Temp0]  =   S[kJb]; 
		S[kS_Temp0] *=   C[kMu_c_l]; 
 
		T4[kT4_Temp0].SumOf ( T4[kAAe_2bar], T4[kEE_2bar] );
		T4[kT4_Temp0] *=  S[kS_Temp0]; 
		B[kBa_Zeta_curl_sE_T].MultAB ( T4[kT4_Temp0], 	B[kB_1hat] );  // Needs neg sign

		T4[kT4_Temp0].DiffOf ( T4[kAAf_2bar], T4[kEE_2bar] );
		T4[kT4_Temp0] *=  S[kS_Temp0]; 
		B[kBb_Zeta_curl_sE_T].MultAB ( T4[kT4_Temp0], 	B[kB_1hat] ); 

	}	
	else {
		B[kBa_Zeta_curl_sE_T] = 0.0; 
		B[kBb_Zeta_curl_sE_T] = 0.0; 
	}

	//-- Summation of all contributing B terms

	B[kBa_Z]  = B[kBa_Zeta_Fbi];
	B[kBa_Z] -= B[kBa_Zeta_Jb];
	B[kBa_Z] -= B[kBa_Zeta_curl_sE_T];

	B[kBb_Z]  = B[kBb_Zeta_Jb];
	B[kBb_Z] -= B[kBb_Zeta_Fbi];
	B[kBb_Z] += B[kBb_Zeta_curl_sE_T];

}

