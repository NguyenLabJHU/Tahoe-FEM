// $Id: MFGP_Bal_EqT.cpp,v 1.8 2004-08-20 05:02:05 raregue Exp $
#include "MFGP_Bal_EqT.h" 

using namespace Tahoe;

/* constructor */
/*
MFGP_Bal_EqT::MFGP_Bal_EqT(void):
	fData_Pro_Displ(NULL),
	fData_Pro_Plast(NULL)
	*/
MFGP_Bal_EqT::MFGP_Bal_EqT(void)
{
	
}

/* destructor */
MFGP_Bal_EqT::~MFGP_Bal_EqT(void) { }
/*
{
	delete fData_Pro_Displ;
	delete fData_Pro_Plast;
}
*/



void MFGP_Bal_EqT::Initialize (int &curr_ip, D3MeshFreeShapeFunctionT *Shapes_displ, D3MeshFreeShapeFunctionT *Shapes_plast, 
							GRAD_MRSSKStV *GRAD_MR_Plast_Mat,
							int &fTime_Step, double fdelta_t) 
{
	n_en_displ  = Shapes_displ->Derivatives_U(curr_ip).MajorDim(); //??
	n_en_plast  = Shapes_plast->Derivatives_U(curr_ip).MajorDim();
	n_sd    	= Shapes_displ->Derivatives_U(curr_ip).MinorDim();	
	n_sd_x_n_sd = n_sd * n_sd;
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	/*
	fData_Pro_Displ = new MFGP_MFA_Data_Processor_DisplT;	
	fData_Pro_Displ->Initialize ( Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) ); //??
	fData_Pro_Plast = new MFGP_MFA_Data_Processor_PlastT;
	fData_Pro_Plast->Initialize ( Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) ); //??
	*/
	Data_Pro_Displ.Initialize ( Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) ); //??
	Data_Pro_Plast.Initialize ( Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) ); //??

	stress = GRAD_MR_Plast_Mat->s_ij();
	moduli = GRAD_MR_Plast_Mat->c_ijkl();
	
	Form_C_List		( GRAD_MR_Plast_Mat );
	Form_B_List		(  );
	//Form_V_S_List	( np1 );
	//Form_VB_List	(  );
}

//---------------------------------------------------------------------

void MFGP_Bal_EqT::Form_LHS_Klambda_Ku(dMatrixT &Klambda, dMatrixT &Ku )  
{
	n_rows_matrix = Klambda.Rows();
	n_cols_matrix = Klambda.Cols();
	dMatrixT Klambdatemp (n_rows_matrix,n_cols_matrix);
	Klambda.MultATBC (B1_d, Culam1, phi_lam  );
	Klambdatemp.MultATBC ( B1_d, Culam2, B4_lam );
	Klambda += Klambdatemp;	 	
		
	n_rows_matrix = Ku.Rows();
	n_cols_matrix = Ku.Cols();
	dMatrixT Kutemp (n_rows_matrix,n_cols_matrix);
	Ku.MultATBC ( B1_d, Cuu1, B1_d );
	Kutemp.MultATBC ( B1_d, Cuu2, B3_d );
	Ku += Kutemp;
}


void MFGP_Bal_EqT::Form_RHS_F_int ( dArrayT &F_int ) 
{
	B1_d.MultTx( stress, F_int );
}


void MFGP_Bal_EqT::Form_B_List (void)
{
		B1_d.Dimension 	( n_sd, n_en_displ);
		B3_d.Dimension 	( n_sd, n_en_displ);
		int dum=1;
		phi_lam.Dimension ( dum, n_en_plast);
		B4_lam.Dimension ( dum, n_sd_x_n_en_plast);
		//NTS: check the allocation of phi_lam and B4_lam
		//B_gradu.Dimension ( dum, n_sd); //B3??	
		
		/*
		fData_Pro_Displ->Set_B1 (B1_d);
		fData_Pro_Displ->Set_B3 (B3_d);
		fData_Pro_Plast->Set_phi (phi_lam);
 		fData_Pro_Plast->Set_B4 (B4_lam); //B4 is scalar  ????
 		*/
 		Data_Pro_Displ.Set_B1 (B1_d);
		Data_Pro_Displ.Set_B3 (B3_d);
		Data_Pro_Plast.Set_phi (phi_lam);
 		Data_Pro_Plast.Set_B4 (B4_lam); //B4 is scalar  ????
}


void MFGP_Bal_EqT::Form_C_List (GRAD_MRSSKStV *GRAD_MR_Plast)
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



