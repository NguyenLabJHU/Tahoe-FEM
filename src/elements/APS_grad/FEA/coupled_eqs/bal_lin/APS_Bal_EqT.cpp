// $Id: APS_Bal_EqT.cpp,v 1.16 2003-10-09 21:46:52 raregue Exp $
#include "APS_Bal_EqT.h" 

using namespace Tahoe;

APS_Bal_EqT::APS_Bal_EqT ( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl, APS_MaterialT *APS_Matl,
								APS_VariableT &np1, APS_VariableT &n, 
								int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes, Shear_Matl, APS_Matl, np1, n, Integration_Scheme);
}

/* destructor */
//APS_Bal_EqT::~APS_Bal_EqT(void) { }


//---------------------------------------------------------------------

void APS_Bal_EqT::Construct ( FEA_ShapeFunctionT &Shapes, APS_MaterialT *Shear_Matl,
							APS_MaterialT *APS_Matl, APS_VariableT &np1, APS_VariableT &n, 
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
	
	//tmp
	n_en_surf=2;

	delta_t = fdelta_t;
	
	Data_Pro.Construct ( Shapes.dNdx	);
	Data_Pro.Insert_N  ( Shapes.N );
	Integral.Construct ( Shapes.j, Shapes.W ); 
	
	Form_C_List		( Shear_Matl, APS_Matl );
	Form_B_List		(  );
	Form_V_S_List	( np1 );
	Form_VB_List	(  );
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Keps_Kd	( dMatrixT &Keps, dMatrixT &Kd )  
{
		Keps 	= Integral.of( B_d[kB], C[kMu], B_eps[kBgamma] );  
		Keps 	*= -1.0;
	 	Kd  	= Integral.of( B_d[kB], C[kMu], B_d[kB] );  	
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int, APS_VariableT &npt ) 
{
		//V[kgrad_u] = npt.Get ( APS::kgrad_u );
		B_gradu[kgrad_u] = npt.Get ( APS::kgrad_u );
		V[kV_Temp2](0)=B_gradu[kgrad_u](0,0);
		V[kV_Temp2](1)=B_gradu[kgrad_u](0,1);
		V[kgammap] = npt.Get ( APS::kgammap );
		//F_int = Integral.of( B[kB], C[kMu], V[kgrad_u] ); 
		F_int = Integral.of( B_d[kB], C[kMu], V[kV_Temp2] ); 
		F_int -= Integral.of( B_d[kB], C[kMu], V[kgammap] );
}


//---------------------------------------------------------------------

void APS_Bal_EqT::Form_LHS_Kd_Surf	( dMatrixT &Kd, FEA_SurfShapeFunctionT &SurfShapes, iArrayT& face_equations )  
{	
		Data_Pro_Surf.Construct ( SurfShapes.dNdx	);
		Data_Pro_Surf.Insert_N_surf  ( SurfShapes.N );
		SurfIntegral.Construct ( SurfShapes.j, SurfShapes.W );
		
		Data_Pro_Surf.APS_B_surf(B_d_surf[kB_surf]);
		Data_Pro_Surf.APS_N(VB_d[kN]);

		V[knueps] = SurfShapes.normal;
		//V[knueps](1) = SurfShapes.normal[1];

 		V[knueps].Dot( B_d_surf[kB_surf], VB_d[knuB] ); 
 		
 		dMatrixT Kd_surfdof(n_en_surf), Kd_eldof(n_en);
		Kd_surfdof	= SurfIntegral.of( VB_d[kN], C[kMu], VB_d[knuB] );
		Kd_surfdof	*= -1.0;
		
		//put Kd_tmp in proper eqnos of Kd
		if (face_equations.Length() == 2)
		{
		for (int i=0; i<face_equations.Length(); i++)
			for (int j=0; j<face_equations.Length(); j++)
				Kd_eldof[face_equations[i]][face_equations[j]] = Kd_surfdof[i][j];				
		}
		Kd += Kd_eldof;
}

//---------------------------------------------------------------------

void APS_Bal_EqT::Form_RHS_F_int_Surf ( dArrayT &F_int, APS_VariableT &npt, double &wght, iArrayT& face_equations  ) 
{
		V[keps](0) = C[km1];
		V[keps](0) *= wght;
		V[keps](1) = C[km2];
		V[keps](1) *= wght;
		
		B_gradu[kgrad_u_surf] = npt.Get ( APS::kgrad_u_surf );
		V[kV_Temp2](0)=B_gradu[kgrad_u_surf](0,0);
		V[kV_Temp2](1)=B_gradu[kgrad_u_surf](0,1);
		V[knueps].Dot( V[kV_Temp2], S[knuepsgradu] );
		V[knueps].Dot( V[keps], S[knuepseps] );

		dArrayT F_int_surfdof(n_en_surf), F_int_eldof(n_en);
		F_int_surfdof = SurfIntegral.of( VB_d[kN], C[kMu], S[knuepsgradu] ); 
		F_int_surfdof *= -1.0;
		F_int_surfdof += SurfIntegral.of( VB_d[kN], C[kMu], S[knuepseps] );
		
		//put components of F_int_tmp in proper eqnos of F_int
		if (face_equations.Length() == 2)
		{
		for (int i=0; i<face_equations.Length(); i++)
			F_int_eldof[face_equations[i]] = F_int_surfdof[i];
		}
		F_int += F_int_eldof;
}



//=== Private =========================================================
	             				
void APS_Bal_EqT::Form_B_List (void)
{
		B_d.Construct 	( kNUM_B_d_TERMS, n_ip, n_sd, n_en);
		B_d_surf.Construct 	( kNUM_B_d_surf_TERMS, n_ip, n_sd, n_en_surf);
		B_eps.Construct ( kNUM_B_eps_TERMS, n_ip, n_sd, n_sd_x_n_en);
		int dum=1;
		//B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, n_sd, dum);
		B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, dum, n_sd);		
		
		Data_Pro.APS_B(B_d[kB]);
 		Data_Pro.APS_Ngamma(B_eps[kBgamma]);
}

		
void APS_Bal_EqT::Form_VB_List (void)
{
		VB_d.Construct 	( kNUM_VB_d_TERMS, 	n_ip, n_en_surf	);
		VB_eps.Construct ( kNUM_VB_eps_TERMS, 	n_ip, n_sd_x_n_en 	);				
}


void APS_Bal_EqT::Form_V_S_List (APS_VariableT &npt)
{
		S.Construct 	( kNUM_S_TERMS, 	n_ip 		);
		V.Construct 	( kNUM_V_TERMS, 	n_ip, n_sd 	);
		int dum=1;
		VS.Construct 	( kNUM_VS_TERMS, 	n_ip, dum 	);
}


void APS_Bal_EqT::Form_C_List (APS_MaterialT *Shear_Matl, APS_MaterialT *APS_Matl)
{
		C.Dimension 	( kNUM_C_TERMS );
		
		C[kMu] 	= Shear_Matl -> Retrieve ( Shear_MatlT::kMu );
		C[km1]  = APS_Matl -> Retrieve ( APS_MatlT::km1	);
		C[km2]  = APS_Matl -> Retrieve ( APS_MatlT::km2	);
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


