// $Id: MFGP_Bal_EqT.cpp,v 1.1 2004-06-18 19:20:15 kyonten Exp $
#include "MFGP_Bal_EqT.h" 

using namespace Tahoe;

MFGP_Bal_EqT::MFGP_Bal_EqT (ShapeFunctionT &Shapes_displ, ShapeFunctionT &Shapes_plast, 
							MFGP_MaterialT *GRAD_MR_Plast_Mat, 
							MFGP_VariableT &np1, MFGP_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Construct (Shapes_displ, Shapes_plast, GRAD_MR_Plast_Mat, np1, n, Integration_Scheme);
}

/* destructor */
//MFGP_Bal_EqT::~MFGP_Bal_EqT(void) { }


//---------------------------------------------------------------------

void MFGP_Bal_EqT::Construct (FEA_ShapeFunctionT &Shapes_displ, FEA_ShapeFunctionT &Shapes_plast, 
							MFGP_MaterialT *GRAD_MR_Plast_Mat,
							MFGP_VariableT &np1, MFGP_VariableT &n, 
							int &fTime_Step, double fdelta_t, int Integration_Scheme) 
{
	Time_Integration_Scheme = Integration_Scheme;

	n_ip = np1.fVars_vector[0].IPs(); 
	/*n_rows_vector = np1.fVars_vector[0].Rows(); 
	n_rows_matrix = np1.fVars_matrix[0].Rows(); 
	n_cols_matrix = np1.fVars_matrix[0].Cols(); 
	*/
	//#pragma message("MFGP_Bal_EqT::Construct: FEA_dVectorT has no Cols() function")

	n_en_displ    	= Shapes_displ.Dphi.Cols();
	n_en_plast    	= Shapes_plast.Dphi.Cols();
	n_sd    	= Shapes_displ.Dphi.Rows();
	//n_sd 		= n_rows_vector;	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	Data_Pro_Displ.Construct ( Shapes_displ.Dphi	);
	Data_Pro_Displ.Insert_N  ( Shapes_displ.phi );
	Integral.Construct (Shapes_displ.W ); 
	Data_Pro_Plast.Construct ( Shapes_plast.Dphi	);
	Data_Pro_Plast.Insert_N  ( Shapes_plast.phi );
	
	Form_C_List		( GRAD_MR_Plast_Matl );
	Form_B_List		(  );
	//Form_V_S_List	( np1 );
	//Form_VB_List	(  );
}

//---------------------------------------------------------------------

void MFGP_Bal_EqT::Form_LHS_Klambda_Ku(dMatrixT &Klambda, dMatrixT &Ku )  
{
		Klambda = Integral.of( B1_d[kB1], C_ulam1[kModuli], phi_lam[kphi] )
		Klambda += Integral.of( B1_d[kB1], C_ulam2[kModuli], B4_lam[kB4] );
		  
	 	Ku  = Integral.of( B1_d[kB1], C_uu1[kModuli], B1_d[kB1] )
	 	Ku  += Integral.of( B1_d[kB1], C_uu2[kModuli], B3_d[kB3] );  	
}

//---------------------------------------------------------------------
/*
void MFGP_Bal_EqT::Form_LHS_Ku_Klambda(dMatrixT &Klamu, dMatrixT &Klamlam )  
{
		Ku = Integral.of( phi_lam[kphi], C_lamu1[kModuli], B1_d[kB1] )
		Ku += Integral.of( phi_lam[kphi], C_lamu2[kModuli], B3_d[kB3] );  
		
	 	Klambda  = Integral.of( phi_lam[kphi], C_lamlam1[kModuli], phi_lam[kphi] )
	 	Klambda  += Integral.of( phi_lam[kphi], C_lamlam2[kModuli], B4_lam[kB4] );  	
}
*/
//---------------------------------------------------------------------

void MFGP_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int, MFGP_VariableT &npt ) 
{
		F_int = Integral.of( B1_d[kB1], fCurrMaterial->s_ij()); 
		//updated stress passed from the constitutive model
}



//=== Private =========================================================
	             				
void MFGP_Bal_EqT::Form_B_List (void)
{
		B1_d.Construct 	( kNUM_B1_d_TERMS, n_ip, n_sd, n_en_displ);
		B3_d.Construct 	( kNUM_B3_d_TERMS, n_ip, n_sd, n_en_displ);
		phi_lam.Construct ( kNUM_phi_lam_TERMS, n_ip, n_sd, n_en_plast);
		B4_lam.Construct ( kNUM_B4_lam_TERMS, n_ip, n_sd, n_sd_x_n_en_plast);
		int dum=1;
		//NTS: check the allocation of phi_lam and B4_lam
		B_gradu.Construct ( kNUM_B_gradu_TERMS, n_ip, dum, n_sd); //B3??	
		
		Data_Pro_Displ.MFGP_B1(B1_d[kB1]);
		Data_Pro_Displ.MFGP_B3(B3_d[kB3]);
		Data_Pro_Plast.MFGP_phi(phi_lam[kN_lam]);
 		Data_Pro_Plast.MFGP_B4(B4_lam[kB4]); //B4 is scalar
}



void MFGP_Bal_EqT::Form_C_List (MFGP_MaterialT *GRAD_MR_Plast)
{
		//C.Dimension( kNUM_C_TERMS );
		int i,j;
		Cuu1.Dimension(6,6); Cuu2.Dimension(6,6);
		Culam1.Dimension(6,1); Culam2.Dimension(6,1);
		Clamu1.Dimension(1,6); Clamu2.Dimension(1,6);
		//Clamlam1.Dimension(1); Clamlam2.Dimension(1);
		
		C[kMu] 	= Shear_Matl -> Retrieve ( Shear_MatlT::kMu );
		C[km1_x]  = APS_Matl -> Retrieve ( APS_MatlT::km1_x	);
		C[km1_y]  = APS_Matl -> Retrieve ( APS_MatlT::km1_y	);
		C[km2_x]  = APS_Matl -> Retrieve ( APS_MatlT::km2_x	);
		C[km2_y]  = APS_Matl -> Retrieve ( APS_MatlT::km2_y	);
		//retrive Cgep/fmoduli, and form the 8 C matrices
		//pass the fModuli from the constitutive model
		for (i=0; i<=7; ++i)
		{
		 for (j=0; j<=14; ++j)
		 {
		  if (i<=5 & j<=5)
		  {
		   Cuu1(i,j)=fModuli(i,j);
		  }
		  if (i<=5 & j>5 & j<=11)
		  {
		   Cuu2(i,j-6)=fModuli(i,j);
		  }
		  if (i<=5 & j>11 & j<=12)
		  {
		   Culam1(i,j-12)=fModuli(i,j);
		  }
		  if (i<=5 & j>12)
		  {
		   Culam2(i,j-13)=fModuli(i,j);
		  }
		  if (i>5 & j<=5)
		  {
		   Clamu1(i-6,j)=fModuli(i,j);
		  }
		  if (i>5 & j>5 & j<=11)
		  {
		   Clamu2(i-6,j-6)=fModuli(i,j);
		  }
		  if (i>5 & j>11 & j<=12)
		  {
		   Clamlam1(i-6,j-12)=fModuli(i,j);
		  }
		  if (i>5 & j>12)
		  {
		   Clamlam2(i-6,j-13)=fModuli(i,j);
		  }
		 }
		}
}



