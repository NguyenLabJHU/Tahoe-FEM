//DEVELOPMENT


#include "FEA.h" 
#include "VMS.h" // <-- Switch name to VMF later 

using namespace Tahoe;

VMF_Virtual_Work_EqT::VMF_Virtual_Work_EqT ( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl,VMS_VariableT &np1,VMS_VariableT &n, 
								int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,Iso_Matl,np1,n,Integration_Scheme);
}

/* destructor */
//VMF_Virtual_Work_EqT::~VMF_Virtual_Work_EqT(void)
//{
//}

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Construct ( 	FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl,VMS_VariableT &np1,VMS_VariableT &n, 
			int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	n_ip 		= np1.fVars[0].IPs(); 
	n_rows	= np1.fVars[0].Rows(); 
	n_cols	= np1.fVars[0].Cols();
	n_en    = Shapes.dNdx.Cols();
  n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	if ( fTime_Step != time_step) { 	// New time step
		// Can put flow rule here: Eventually do this with all the kVar_n and make method Next_Step()
		time_step = fTime_Step;
	}

	delta_t = fdelta_t;
	
	Data_Pro.Construct ( Shapes.dNdx );

	Form_C_List			(	Iso_Matl );
  Form_A_S_Lists 	(	np1,n		 );
  Form_B_List			(	);

	//A.Print("A Matricies");
	//T4.Print("T4 Matrix");
	//S.Print("S Scalar");
	//B.Print("B Matricies");

	Integral.Construct ( Shapes.j, Shapes.W ); 

}

//---------------------------------------------------------------------

/** RECALL CLASSIC NEWTON-RAPHSON :
 
    K.delta_d = -[k.d - f]  : Where both K and k are element matricies.
    And where K is the tangent, drive it to zero for roots (RHS and delta_d will also vanish) 
		Note, k.d are the internal forces (F_int) and f are the external forces (F_ext).

	  MULTI-FIELD NEWTON-RAPHSON :	

    Ka.delta_da + Kb.delat_db = -[F_int - f]  : Where Ka,Kb, and k are element matricies.
    Ka and Kb are tangents, drive them to zero in a staggered scheme (RHS, delta_da, and delta_db
	 	will also vanish).  Note, dislocation glide (da) contributes to stress in a round-a-bout way.
		While it is true that the intermediate configuration is stress free, second PK S = S(Fb) 
		and Fb := dx/dX_bar. Recall X_bar = X + ua.  So changes in ua will ultimately affect S and
		sigma sinse sigma = 1/j FbSFb^T.  */	

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Form_LHS_Ka_Kb	( dMatrixT &Ka, dMatrixT &Kb )  // Untested
{
	/* Term I. 		*/		Ka 	= Integral.of( 	B[kB_1hat], B[kBI_tau_3hat] 									);  	
	/* Term IIb. 	*/	 	Kb  = Integral.of( 	B[kB_1hat], B[kBbII_2hat] 										);  	
	/* Term IIa. 	*/	 	Ka -= Integral.of( 	B[kB_1hat], B[kBaII_3hat] 										);  	
	/* Term IIIb.	*/	 	Kb += Integral.of( 	B[kB_1hat], T4[kcc_b], 				B[kB_1hat] 			);  	
	/* Term IIIa.	*/	 	Ka -= Integral.of( 	B[kB_1hat], T4[kcc_b], 				B[kBaIII_2bar] 	);  	
}

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Form_RHS_F_int ( dArrayT &F_int ) // Untested
{
	FEA_dVectorT sigma_vec	( n_ip, n_sd_x_n_sd 		); // <-- Dimensionality problem
	Data_Pro.Reduce_Order		(	A[kSigma], sigma_vec 	); 

	F_int = Integral.of			( B[kB_1hat], sigma_vec	);  // <-- sigma_vec must be dim n_sd_x_n_en	
}

//=== Private =========================================================

//##################################################################################
//################ B_TERMS #########################################################
//##################################################################################

void VMF_Virtual_Work_EqT::Form_B_List (void)
{
		B.Construct ( kNUM_B_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_en );	

	 	Data_Pro.grad_u       		( 														B[kB_1hat]  			); 
	 	Data_Pro.A_grad_u_T_B 		( 	A[kSigma], 		A[kFbT],  	B[kBI_tau_3hat] 	);
	 	Data_Pro.grad_u_A			 		( 	A[kSigma], 								B[kBbII_2hat] 		);
	 	Data_Pro.A_grad_u_B			 	( 	A[kF_sharp],	A[kSigma], 	B[kBaII_3hat] 		);
	 	Data_Pro.A_grad_u				 	( 	A[kF_sharp],						 	B[kBaIII_2bar] 		);
}
						
//##################################################################################
//################ A_TERMS #########################################################
//##################################################################################

//NOTE: np1 := "n+1" time step; n := "n" time step; npt := "n+theta" time step
//      *** No subscript implies n+theta time step ***

void VMF_Virtual_Work_EqT::Form_A_S_Lists (VMS_VariableT &npt,VMS_VariableT &n,int Integration_Scheme)
{

#if 0
	//---- Developer cheat: put npt in function door in-lieu-of np1 or n for speed 
	VMS_VariableT npt(	n.Get(VMS::kGRAD_ua), n.Get(VMS::kGRAD_ub)	);  

	if 			(	Integration_Scheme == FEA::kForward_Euler		)		npt = n;
	else if (	Integration_Scheme == FEA::kBackward_Euler	)		npt = np1;
	else if (	Integration_Scheme == FEA::kCrank_Nicholson	) { npt.SumOf(np1,n); npt *= 0.5; }
	else 	cout << " ...ERROR >> VMF_Virtual_Work_EqT::Form_A_List() : Bad theta value for time stepping \n";
#endif

	A.Construct ( kNUM_A_TERMS, n_ip, n_sd, n_sd);
	S.Construct ( kNUM_S_TERMS, n_ip);

  A[kF]   		= npt.Get (	VMS::kF					);	 // NOTE: kF != VMS::kF 	
  A[kFb]  		= npt.Get (	VMS::kFb				);	 	
	A[kgrad_ub]	= npt.Get (	VMS::kgrad_ub		);
	S[kJb].Determinant 		( A[kFb] 					);

	//--- For CCba GammaY
	A[kFb_n]  = n.Get			(	VMS::kFb );
	A[kCb_n].MultATB			( A[kFb_n], 	A[kFb_n] );

	A[kF].Determinant	 	  ( S[kJ] 								);  
	A[kFbT].Transpose  		( A[kFb] 								);
	A[kF_sharp].MultAB  	( A[kgrad_ub], 	A[kFb] 	);

	//----- Calculate Eb 
	
	A[kCb].MultATB   			( A[kFb],  			A[kFb] 	);					
	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity		( -1.0 ); 
	A[kEb] *= 0.5; 

	//----- Calculate stresses S and Sigma

	T4.Construct ( kNUM_T4_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_sd );  // Move to Initialize() later

 	if 	(bControl_Eb)  {
		Form_CC ( );  // CC becomes either classic CC or CCvmst, Lame_ba formed if needed
  	Data_Pro.C_IJKL_E_KL	(  	S[kLamda_ba], 	S[kMu_ba], 	A[kEb], 	A[kS] 		); // UT
		A[kSigma].MultABCT 		(  	A[kFb], 				A[kS], 			A[kFb] 							);
		A[kSigma] /= S[kJ]; // not kJb !!
  	Data_Pro.c_ijkl				(		S[kJ], 					A[kFb], 		T4[kCC], 	T4[kcc_b]	); // This Works (use for non-sym CC) : Slow
	}
	else {
  	Data_Pro.C_IJKL 			(		C[kLamda],  C[kMu],		T4[kCC]	  			);  
  	Data_Pro.C_IJKL_E_KL	(  	C[kLamda], 	C[kMu], 	A[kEb],  A[kS] 	); // UT
		A[kSigma].MultABCT 		(  	A[kFb], 		A[kS], 		A[kFb] 					);
		A[kSigma] /= S[kJb];  // kJb can be kJ (if Plastic Volume is conserved !!
  	Data_Pro.c_ijkl				(		C[kLamda], 	C[kMu], S[kJ], A[kFb], T4[kcc_b]	); // This Works (original) 
	}

	//----- 4th Order VMF Finite Strain Elasticity Tensor
  	
	//Data_Pro.c_ijkl				(	S[kLamda_ba], S[kMu_ba], S[kJ], A[kFb], T4[kcc_b]	); // Was original method (uses symmetry)
  // Data_Pro.c_ijkl				(	C[kLamda], C[kMu], S[kJ], A[kFb], T4[kcc_b]	); // This Works (original) 
	// The following can be used for testing but don't work with a straight substitution
  //Data_Pro.c_ijkl_Alt		(	lamda,mu, S[kJ], A[kFb], T4[kT4_Temp0]	); T4[] is smaller due to use of symmetry (test only)
  //Data_Pro.c_ijkl				(	T4[kCC], S[kJ], A[kFb], T4[kcc_b]	);	// cc is smaller, expand before implementation 
 	

}

//##################################################################################
//################ C_TERMS #########################################################
//##################################################################################

void VMF_Virtual_Work_EqT::Form_C_List (VMF_MaterialT *Iso_Matl)
{
	 C.Dimension 	(	kNUM_C_TERMS );

	 C[kLamda]	  = Iso_Matl -> Retrieve ( Iso_MatlT::kLamda 		);
	 C[kMu] 			= Iso_Matl -> Retrieve ( Iso_MatlT::kMu 			);

	 //-- Cba Constants
	 C[kLamda1]  	= Iso_Matl -> Retrieve ( Iso_MatlT::kLamda1 	);
	 C[kMu1]  		= Iso_Matl -> Retrieve ( Iso_MatlT::kMu1 			);
	 C[kLamda2]  	= Iso_Matl -> Retrieve ( Iso_MatlT::kLamda2 	);
	 C[kMu2]  		= Iso_Matl -> Retrieve ( Iso_MatlT::kMu2 			);
	 C[kPi] 			= Iso_Matl -> Retrieve ( Iso_MatlT::kPi 			);
	 C[kRho]  		= Iso_Matl -> Retrieve ( Iso_MatlT::kRho 			);
	 C[kGamma_b]  = Iso_Matl -> Retrieve ( Iso_MatlT::kGamma_b 	);
	 C[kAlphaY]  	= Iso_Matl -> Retrieve ( Iso_MatlT::kAlphaY 	);
	 C[k2Gamma_b]  = 2.0 * C[kGamma_b]; 
	 C[kAlphaY_plus_1]  =  C[kAlphaY] + 1.0; 
 
}


//################################## CCba ################################################

//---------------------------------------------------------------------------------------
//NOTE: the only place CC is used is in forming the LHS tangent, do it makes sense
//to pack CCvmst into CC.  Recall that on the RHS, S is formed by using Lame Constants.

void VMF_Virtual_Work_EqT::Form_CC (void) //-- Puts either regular CC or CCvmst in CC 
{

		CCba_Scalars ( );  // Also creates Eb_dot and N_Eb_dot, Mag(Eb_dot)

		//cout << "Lamda1, Mu1 = "<<C[kLamda1]<<", "<<C[kMu1]<<endl;

		//-- Build CC1 & CC2 
  	Data_Pro.C_IJKL (	C[kLamda1],  C[kMu1],		T4[kCC1]	); 
  	Data_Pro.C_IJKL (	C[kLamda2],  C[kMu2],		T4[kCC2]	);  

		//-- Build CC_ba
		T4[kCC_ba]  = T4[kCC1];
		T4[kCC_ba] *= S[kPsi1];
		T4[kCC_ba] += T4[kCC2]; 
	
		//-- Build CC_gamma	
		Data_Pro.A_o_B 					( A[kEb], 			A[kN_Eb_dot],	T4[kEEn] 			);  
		T4[kCC_gamma].MultAB		( T4[kCC1],			T4[kEEn]										);
		T4[kCC_gamma] *= S[kPsi]; 
		T4[kCC_gamma] /= ( delta_t * C[kGamma_b] ); 
		T4[kCC_gamma] *= C[kAlphaY_plus_1]; 

		//-- Build CC_psi	
		A[kA_Temp0] = A[kEb];
		Data_Pro.A_o_B 					( A[kA_Temp0], 	A[kEb],				T4[kEb_o_Eb] 	);  
		T4[kCC_psi].MultAB			( T4[kCC1],			T4[kEb_o_Eb]								);
		T4[kCC_psi] *= S[kGammaY_bar_Eta1_plus_Eta2];

		//-- Build CCvmst (pack it in CC, go on as usual with del(S) term)
		T4[kCC]  = T4[kCC_psi];
		T4[kCC] += T4[kCC_ba];
		T4[kCC] += T4[kCC_gamma]; // <-- Comment/Uncomment (Turn term off/on) for no rate dependency

}


//################################## CCba ################################################

void VMF_Virtual_Work_EqT::CCba_Scalars( void )  
{
  A[kEb].Mag_and_Dir ( S[kMag_Eb], A[kNeb] ); 

	if ( S[kMag_Eb] != 0.0 ) {

  	S[kRho_Mag_Eb]  = S[kMag_Eb]; 
  	S[kRho_Mag_Eb] *= C[kRho]; 

		S[kPi_Tanh].Tanh( S[kRho_Mag_Eb] ); 
		S[kPi_Tanh] *= C[kPi];
		S[kPsi]  	 	 = S[kPi_Tanh]; 
		S[kPsi] 		/= S[kMag_Eb]; 

		S[kEta1]  = S[kPsi];			// Eta1 = Psi/(Mag_Eb^2)
		S[kEta1] /= S[kMag_Eb]; 
		S[kEta1] /= S[kMag_Eb]; 

		S[kEta2].Sech( S[kRho_Mag_Eb] ); 
		S[kEta2].Squared ( );
		S[kEta2] *= C[kPi]; 
		S[kEta2] /= S[kMag_Eb]; 
		S[kEta2] /= S[kMag_Eb]; 

		S[kEta1_plus_Eta2]  = S[kEta1];
		S[kEta1_plus_Eta2] += S[kEta2];

		S[kPsi1] 		= S[kPsi];
		S[kGammaY_bar_Eta1_plus_Eta2]  = S[kEta1_plus_Eta2];

		//------------- Add rate dependency to S_tilde (Modify Psi1 with rate term GammaY)
#if 1 
		//--- Get Eb_dot
		
		A[kEb_n]  = A[kCb_n]; 
		A[kEb_n].PlusIdentity (-1.0); 
		A[kEb_n] *= 0.5; 

		A[kEb_dot].DiffOf ( A[kEb], A[kEb_n] );
		A[kEb_dot] /= delta_t; 
  	A[kEb_dot].Mag_and_Dir 		( S[kMag_Eb_dot], 	A[kN_Eb_dot] );

		//--- Finish Psi1
	
		S[kGammaY]  	  = S[kMag_Eb_dot];	
		S[kGammaY] 		 /= C[kGamma_b];	
		S[kGammaY_bar]  = S[kGammaY];	
		S[kGammaY_bar] *= C[kAlphaY_plus_1];	
		S[kGammaY_bar] -= C[kAlphaY];	

		//-- Comment/Uncomment these for rate dependency
		S[kPsi1]   *= S[kGammaY_bar];
		S[kGammaY_bar_Eta1_plus_Eta2]  *= S[kGammaY_bar];
#endif 
		//-------------- End Modification
		
		//-- Build Lame_ba Constants
	  S[kLamda_ba]  = 	S[kPsi1];
	  S[kLamda_ba] *= 	C[kLamda1];
	  S[kLamda_ba] += 	C[kLamda2];
	  S[kMu_ba]  		= 	S[kPsi1];
	  S[kMu_ba] 	 *= 	C[kMu1];
	  S[kMu_ba] 	 += 	C[kMu2];

	}
	else { // So it doesn't blow up
		S[kPsi1] = 0.0;
		S[kEta1_plus_Eta2]  = 0.0;
		S[kGammaY_bar_Eta1_plus_Eta2]  = 0.0;
	}


}

//##################################################################################

void VMF_Virtual_Work_EqT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "Sigma" )
		tensor = A[kSigma];
	else if ( Name == "S" )
		tensor = A[kS];
	else if ( Name == "Eb" )
		tensor = A[kEb];
	else if ( Name == "Eb_tilde" )
		tensor = A[kEb_tilde];
	else if ( Name == "F" )
		tensor = A[kF];
	else if ( Name == "Fb" )
		tensor = A[kFb];
	else if ( Name == "Cb" )
		tensor = A[kCb];
	else if ( Name == "grad_ub" )
		tensor = A[kgrad_ub];
	else
		cout << " ...ERROR: VMF_Virtual_Work_EqT::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}

//##################################################################################

void VMF_Virtual_Work_EqT::Get ( StringT &Name, FEA_dScalarT &scalar )
{
	if ( Name == "J" )
		scalar = S[kJ];
	else if ( Name == "Jb" )
		scalar = S[kJ];
	else if ( Name == "Rho_Mag_Eb" )
		scalar = S[kRho_Mag_Eb];
	else
		cout << " ...ERROR: VMF_Virtual_Work_EqT::Get() >> Unknown scalar '"<<Name<<"' requested. \n";
}

