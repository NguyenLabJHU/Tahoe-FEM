// $Id: APS_Bal_EqT.cpp,v 1.12 2003-09-27 00:04:02 raregue Exp $
#include "APS_Bal_EqT.h" 

using namespace Tahoe;

APS_Bal_EqT::APS_Bal_EqT ( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_VariableT &np1, APS_VariableT &n, 
								int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes,Shear_Matl,np1,n,Integration_Scheme);
}

/* destructor */
//APS_Bal_EqT::~APS_Bal_EqT(void) { }


//---------------------------------------------------------------------

void APS_Bal_EqT::Construct ( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_VariableT &np1, APS_VariableT &n, 
			int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	n_ip 		= np1.fVars_vector[0].IPs(); 
	/*n_rows_vector = np1.fVars_vector[0].Rows(); 
	n_rows_matrix = np1.fVars_matrix[0].Rows(); 
	n_cols_matrix = np1.fVars_matrix[0].Cols(); 
	*/
//#pragma message("APS_Bal_EqT::Construct: FEA_dVectorT has no Cols() function")

	n_en    	= Shapes.dNdx.Cols();
	n_sd    	= Shapes.dNdx.Rows();
	//n_sd 		= n_rows_vector;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en = n_sd * n_en;

	if ( fTime_Step != time_step) { 	// New time step
		// Can put flow rule here: Eventually do this with all the kVar_n and make method Next_Step()
		time_step = fTime_Step;
	}

	delta_t = fdelta_t;
	
	Data_Pro.Construct ( Shapes.dNdx 	);
	Data_Pro.Insert_N  ( Shapes.N 		);
	Integral.Construct ( Shapes.j, Shapes.W ); 
	
	Form_C_List		( Shear_Matl );
	Form_B_List		(  );
	Form_V_S_List	( np1 );
	Form_VB_List	(  );


}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd )  // Untested
{
		Keps 	= Integral.of( B_d[kB], C[kMu], B_eps[kBgamma] );  
		Keps 	*= -1.0;
	 	Kd  	= Integral.of( B_d[kB], C[kMu], B_d[kB] );  	
#pragma message("APS_Bal_EqT::Form_LHS_Keps_Kd: this domain over Gamma_eps")
//		Kd		-= Integral.of( VB_d[kN], C[kMu], VB_d[knuB] );
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int, APS_VariableT &npt ) // Untested
{
		//V[kgrad_u] = npt.Get ( APS::kgrad_u );
		B_gradu[kgrad_u] = npt.Get ( APS::kgrad_u );
		V[kV_Temp2](0)=B_gradu[kgrad_u](0,0);
		V[kV_Temp2](1)=B_gradu[kgrad_u](0,1);
		V[kgammap] = npt.Get ( APS::kgammap );
		//F_int = Integral.of( B[kB], C[kMu], V[kgrad_u] ); 
		F_int = Integral.of( B_d[kB], C[kMu], V[kV_Temp2] ); 
		F_int -= Integral.of( B_d[kB], C[kMu], V[kgammap] );
#pragma message("APS_Bal_EqT::Form_RHS_F_int: this domain over Gamma_eps")
/*		F_int -= Integral.of( VB_d[kN], C[kMu], S[knuepsgradu] ); 
		F_int += Integral.of( VB_d[kN], C[kMu], S[knuepseps] );  */
}

//=== Private =========================================================
	             				
void APS_Bal_EqT::Form_B_List (void)
{
		B_d.Construct 	( kNUM_B_d_TERMS, n_ip, n_sd, n_en);
		B_eps.Construct ( kNUM_B_eps_TERMS, n_ip, n_sd, n_sd_x_n_en);
		int dum=1;
		//B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, n_sd, dum);
		B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, dum, n_sd);		
		
		Data_Pro.APS_B(B_d[kB]);
 		Data_Pro.APS_Ngamma(B_eps[kBgamma]);
}

		
void APS_Bal_EqT::Form_VB_List (void)
{

		VB_d.Construct 	( kNUM_VB_d_TERMS, 	n_ip, n_en 	);
		VB_eps.Construct ( kNUM_VB_eps_TERMS, 	n_ip, n_sd_x_n_en 	);	
						
		Data_Pro.APS_N(VB_d[kN]);

 		V[knueps].Dot( B_d[kB], VB_d[knuB] );
 		V[knueps].Dot( B_eps[kBgamma], VB_eps[knuNgam] ); 
}


void APS_Bal_EqT::Form_V_S_List (APS_VariableT &npt)
{
		S.Construct 	( kNUM_S_TERMS, 	n_ip 		);
		V.Construct 	( kNUM_V_TERMS, 	n_ip, n_sd 	);
		int dum=1;
		VS.Construct 	( kNUM_VS_TERMS, 	n_ip, dum 	);
		
		#pragma message("APS_Bal_EqT::Form_V_S_List: V[knueps] and V[keps] must be input for BC")
		V[knueps](0) = 1.0;
		V[knueps](1) = 0.0;
		V[keps](0) = 1.0;
		V[keps](1) = 0.0;
		//V[kgrad_u] = npt.Get ( APS::kgrad_u );
		B_gradu[kgrad_u] = npt.Get ( APS::kgrad_u );
		V[kV_Temp2](0)=B_gradu[kgrad_u](0,0);
		V[kV_Temp2](1)=B_gradu[kgrad_u](0,1);
		//V[knueps].Dot( V[kgrad_u], S[knuepsgradu] );
		//V[knueps].Dot( B_gradu[kgrad_u], VS[kVS_Temp1] );
		V[knueps].Dot( V[kV_Temp2], S[knuepsgradu] );
		//S[knuepsgradu]=VS[kVS_Temp1](0);
		V[knueps].Dot( V[keps], S[knuepseps] );
}


void APS_Bal_EqT::Form_C_List (APS_MaterialT *Shear_Matl)
{
		C.Dimension 	( kNUM_C_TERMS );
		
		C[kMu] 	= Shear_Matl -> Retrieve ( Shear_MatlT::kMu );
}


void APS_Bal_EqT::Get ( StringT &Name, FEA_dMatrixT &tensor )
{
	if ( Name == "grad_u" )
		tensor = B_gradu[kgrad_u];
	else
		cout << " ...ERROR: APS_Bal_EqT::Get() >> Unknown tensor '"<<Name<<"' requested. \n";
}


void APS_Bal_EqT::Get ( StringT &Name, FEA_dVectorT &vector )
{
/*	if ( Name == "grad_u" )
		vector = V[kgrad_u]; */
	if ( Name == "gammap" )
		vector = V[kgammap];
	else
		cout << " ...ERROR: APS_Bal_EqT::Get() >> Unknown vector '"<<Name<<"' requested. \n";
}

void APS_Bal_EqT::Get ( StringT &Name, FEA_dScalarT &scalar )
{
//for now, no scalars to get
/*	if ( Name == "?" )
		scalar = S[?];
	else if ( Name == "?" )
		scalar = S[?];
	else
		cout << " ...ERROR: APS_Bal_EqT::Get() >> Unknown scalar '"<<Name<<"' requested. \n";
		*/
}


