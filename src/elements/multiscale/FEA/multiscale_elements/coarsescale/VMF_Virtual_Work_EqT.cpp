//DEVELOPMENT


#include "FEA.h" 
#include "VMS.h" // <-- Switch name to VMF later 

using namespace Tahoe;

VMF_Virtual_Work_EqT::VMF_Virtual_Work_EqT	( FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl,VMS_VariableT &np1,VMS_VariableT &n, 
											int Integration_Scheme) 
{
	Construct (Shapes,Iso_Matl,np1,n,Integration_Scheme);
}

/* destructor */
//VMF_Virtual_Work_EqT::~VMF_Virtual_Work_EqT(void)
//{
//}

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Construct ( 	FEA_ShapeFunctionT &Shapes,VMF_MaterialT *Iso_Matl,VMS_VariableT &np1,VMS_VariableT &n, 
														int Integration_Scheme) 
{
#pragma unused(Integration_Scheme)

	n_ip 		= np1.fVars[0].IPs(); 
	n_rows	= np1.fVars[0].Rows(); 
	n_cols	= np1.fVars[0].Cols();
	n_en    = Shapes.dNdx.Cols();
  n_sd 		= n_rows;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;
	lamda = Iso_Matl -> Retrieve ( Iso_MatlT::kLamda 	);
	mu 		= Iso_Matl -> Retrieve ( Iso_MatlT::kMu 		);
 
	Data_Pro.Construct ( Shapes.dNdx );

  Form_A_S_Lists 	(np1,n);
  Form_B_List			();

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

void VMF_Virtual_Work_EqT::Form_LHS_Ka_Kb	( dMatrixT &Ka, dMatrixT &Kb, double delta_t)  // Untested
{
#pragma unused(delta_t)

	/* Term I. 		*/		Ka 	= Integral.of( 	B[kB_1hat], B[kBI_tau_3hat] 								);  	
	/* Term IIb. 	*/	 	Kb  = Integral.of( 	B[kB_1hat], B[kBbII_2hat] 									);  	
	/* Term IIa. 	*/	 	Ka -= Integral.of( 	B[kB_1hat], B[kBaII_3hat] 									);  	
	/* Term IIIb.	*/	 	Kb += Integral.of( 	B[kB_1hat], T4[kd_1bar],		B[kB_1hat] 			);  	
	/* Term IIIa.	*/	 	Ka -= Integral.of( 	B[kB_1hat], T4[kd_1bar],		B[kBaIII_2bar] 	);  	
}

//---------------------------------------------------------------------

void VMF_Virtual_Work_EqT::Form_RHS_F_int ( dArrayT &F_int, double delta_t) // Untested
{
#pragma unused(delta_t)

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

	A[kF].Determinant	 	  ( S[kJ] 								);  
	A[kFbT].Transpose  		( A[kFb] 								);
	A[kF_sharp].MultAB  	( A[kgrad_ub], 	A[kFb] 	);

	A[kCb].MultATB   			( A[kFb],  			A[kFb] 	);					
	A[kEb]  = A[kCb]; 
	A[kEb].PlusIdentity(-1.0); 
	A[kEb] *= 0.5; 

	//----- Calculate stresses S and Sigma
	
  Data_Pro.C_IJKL_E_KL	(  lamda, mu, A[kEb], A[kS] 	); // Untested Hooke's Law
	A[kSigma].MultABCT 		(  A[kFb], 		A[kS], 	A[kFb] 	);
	A[kSigma] /= S[kJ]; // not kJb !!

	//----- 4th Order VMF Finite Strain Elasticity Tensor

	T4.Construct ( kNUM_T4_TERMS, n_ip, n_sd_x_n_sd, n_sd_x_n_sd );

  Data_Pro.c_ijkl				(	lamda,mu, S[kJ], A[kFb], T4[kd_1bar]	);
  Data_Pro.c_ijkl_Alt		(	lamda,mu, S[kJ], A[kFb], T4[kT4_Temp0]	);

}


