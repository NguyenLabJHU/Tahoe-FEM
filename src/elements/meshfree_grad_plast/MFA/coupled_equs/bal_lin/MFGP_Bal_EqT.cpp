// $Id: MFGP_Bal_EqT.cpp,v 1.4 2004-08-14 00:05:28 raregue Exp $
#include "MFGP_Bal_EqT.h" 

using namespace Tahoe;

/* constructor */
MFGP_Bal_EqT::MFGP_Bal_EqT ( int &curr_ip, D3MeshFreeShapeFunctionT &Shapes_displ, D3MeshFreeShapeFunctionT &Shapes_plast, 
							GRAD_MRSSKStV &GRAD_MR_Plast_Mat, 
							int &fTime_Step, double fdelta_t) 
{
	//Initialize (curr_ip, Shapes_displ, Shapes_plast, GRAD_MR_Plast_Mat, fTime_Step, fdelta_t);
}

/* destructor */
MFGP_Bal_EqT::~MFGP_Bal_EqT(void) { }

void MFGP_Bal_EqT::Initialize (int &curr_ip, D3MeshFreeShapeFunctionT &Shapes_displ, D3MeshFreeShapeFunctionT &Shapes_plast, 
							GRAD_MRSSKStV &GRAD_MR_Plast_Mat,
							int &fTime_Step, double fdelta_t) 
{
	n_en_displ  = Shapes_displ.Dphi.Cols(); //??
	n_en_plast  = Shapes_plast.Dphi.Cols();
	n_sd    	= Shapes_displ.Dphi.Rows();	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	Data_Pro_Displ.Initialize ( Shapes_displ.Derivatives_U(curr_ip), Shapes_displ.DDDerivatives_U(curr_ip) ); //??
	//Data_Pro_Displ.Initialize ( Shapes_displ.Dphi );
	
	Data_Pro_Plast.Initialize ( Shapes_displ.IPShapeU(curr_ip), Shapes_displ.DDerivatives_U(curr_ip)	); //??
	//Data_Pro_Plast.Initialize ( Shapes_plast.Dphi	);
	
	stress = GRAD_MR_Plast.s_ij();
	moduli = GRAD_MR_Plast.c_ijkl();
	
	Form_C_List		( GRAD_MR_Plast_Mat );
	Form_B_List		(  );
	//Form_V_S_List	( np1 );
	//Form_VB_List	(  );
}

//---------------------------------------------------------------------

void MFGP_Bal_EqT::Form_LHS_Klambda_Ku(dMatrixT &Klambda, dMatrixT &Ku )  
{
	n_rows = Klambda.Rows();
	n_cols = Klambda.Cols();
	dMatrixT Klambdatemp (n_rows,n_cols);
	Klambda.MultATBC (B1_d, C_ulam1, phi_lam  );
	Klambda += Klambdatemp.MultATBC ( B1_d, C_ulam2, B4_lam );	 	
		
	n_rows = Ku.Rows();
	n_cols = Ku.Cols();
	dMatrixT Kutemp (n_rows,n_cols);
	Ku.MultATBC ( B1_d, C_uu1, B1_d );
	Ku += Kutemp.MultATBC ( B1_d, C_uu2, B3_d );
}


void MFGP_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int ) 
{
	F_int.MultATB( B1_d, stress);
}


void MFGP_Bal_EqT::Form_B_List (void)
{
		B1_d.Dimension 	( n_sd, n_en_displ);
		B3_d.Dimension 	( n_sd, n_en_displ);
		int dum=1;
		phi_lam.Dimension ( dum, n_en_plast);
		B4_lam.Dimension ( dum, n_sd_x_n_en_plast);
		//NTS: check the allocation of phi_lam and B4_lam
		B_gradu.Dimension ( dum, n_sd); //B3??	
		
		Data_Pro_Displ.Set_B1 (B1_d);
		Data_Pro_Displ.Set_B3 (B3_d);
		Data_Pro_Plast.Set_phi (phi_lam);
 		Data_Pro_Plast.Set_B4 (B4_lam); //B4 is scalar  ????
}


void MFGP_Bal_EqT::Form_C_List (GRAD_MRSSKStV &GRAD_MR_Plast)
{
		//C.Dimension( kNUM_C_TERMS );
		int i,j;
		Cuu1.Dimension(6,6); Cuu2.Dimension(6,6);
		Culam1.Dimension(6,1); Culam2.Dimension(6,1);
		Clamu1.Dimension(1,6); Clamu2.Dimension(1,6);
		//Clamlam1.Dimension(1); Clamlam2.Dimension(1);
		
		//retrive Cgep/fmoduli, and form the 8 C matrices
		//pass the moduli from the constitutive model
		for (i=0; i<=7; ++i)
		{
		 for (j=0; j<=14; ++j)
		 {
		  if (i<=5 & j<=5)
		  {
		   Cuu1(i,j)=moduli(i,j);
		  }
		  if (i<=5 & j>5 & j<=11)
		  {
		   Cuu2(i,j-6)=moduli(i,j);
		  }
		  if (i<=5 & j>11 & j<=12)
		  {
		   Culam1(i,j-12)=moduli(i,j);
		  }
		  if (i<=5 & j>12)
		  {
		   Culam2(i,j-13)=moduli(i,j);
		  }
		  if (i>5 & j<=5)
		  {
		   Clamu1(i-6,j)=moduli(i,j);
		  }
		  if (i>5 & j>5 & j<=11)
		  {
		   Clamu2(i-6,j-6)=moduli(i,j);
		  }
		  if (i>5 & j>11 & j<=12)
		  {
		   Clamlam1(i-6,j-12)=moduli(i,j);
		  }
		  if (i>5 & j>12)
		  {
		   Clamlam2(i-6,j-13)=moduli(i,j);
		  }
		 }
		}
}



