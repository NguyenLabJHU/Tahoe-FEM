// $Id: MFGP_Bal_EqT.cpp,v 1.12 2004-10-24 17:35:40 kyonten Exp $
#include "MFGP_Bal_EqT.h" 

using namespace Tahoe;

/* constructor */
/*
MFGP_Bal_EqT::MFGP_Bal_EqT(void):
	fData_Pro_Displ(NULL),
	fData_Pro_Plast(NULL)
	*/
MFGP_Bal_EqT::MFGP_Bal_EqT(void) { }


/* destructor */
//MFGP_Bal_EqT::~MFGP_Bal_EqT(void) { }
/*
{
	delete fData_Pro_Displ;
	delete fData_Pro_Plast;
}
*/


void MFGP_Bal_EqT::Initialize(int& curr_ip, D3MeshFreeShapeFunctionT *Shapes_displ, D3MeshFreeShapeFunctionT *Shapes_plast, 
							GRAD_MRSSKStV *GRAD_MR_Plast_Mat,
							int& fTime_Step, double fdelta_t) 
{
	n_en_displ  = Shapes_displ->Derivatives_U(curr_ip).MinorDim(); //dN: [nsd]x[nnd]
	n_en_plast  = Shapes_plast->Derivatives_U(curr_ip).MinorDim();
	n_sd    	= Shapes_displ->Derivatives_U(curr_ip).MajorDim();
	n_str = dSymMatrixT::NumValues(n_sd);	
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
	
	/*
	fData_Pro_Displ = new MFGP_MFA_Data_Processor_DisplT;	
	fData_Pro_Displ->Initialize ( Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) ); //??
	fData_Pro_Plast = new MFGP_MFA_Data_Processor_PlastT;
	fData_Pro_Plast->Initialize ( Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) ); //??
	*/
	Data_Pro_Displ.Initialize(Shapes_displ->Derivatives_U(curr_ip), Shapes_displ->DDDerivatives_U(curr_ip) ); //??
	Data_Pro_Plast.Initialize(Shapes_displ->IPShapeU(curr_ip), Shapes_displ->DDerivatives_U(curr_ip) ); //??

	stress = GRAD_MR_Plast_Mat->s_ij();
	moduli = GRAD_MR_Plast_Mat->c_ijkl();
	
	Form_C_List(GRAD_MR_Plast_Mat);
	Form_B_List(  );
}

//---------------------------------------------------------------------

void MFGP_Bal_EqT::Form_LHS_Klambda_Ku(dMatrixT& Klambda, dMatrixT& Ku)  
{
	n_rows_matrix = Klambda.Rows();
	n_cols_matrix = Klambda.Cols();
	dMatrixT Klambdatemp(n_rows_matrix,n_cols_matrix);
	Klambda.MultATBC(B1, Culam1, psi_lam);
	Klambdatemp.MultATBC(B1, Culam2, B4);
	Klambda += Klambdatemp;	// Klambda: [nsd*nnd]x[nnd] 	
		
	n_rows_matrix = Ku.Rows();
	n_cols_matrix = Ku.Cols();
	dMatrixT Kutemp(n_rows_matrix,n_cols_matrix);
	Ku.MultATBC(B1, Cuu1, B1);
	Kutemp.MultATBC(B1, Cuu2, B3);
	Ku += Kutemp; // Ku: [nsd*nnd]x[nsd*nnd]
}


void MFGP_Bal_EqT::Form_RHS_F_int(dArrayT& F_int) 
{
	
	dMatrixT B1_T(n_sd_x_n_en_displ, n_str);
	dArrayT stress_sym(n_str);
	B1_T.Transpose(B1);
	
	/* take out the symmetric part of the stress */
	for (int ij = 0; ij < n_str; ij++)
	{
		int i, j;
		dSymMatrixT::ExpandIndex(n_sd, ij, i, j);
		stress_sym[ij] = stress(i,j);
	} 
	
	B1_T.MultTx(stress_sym, F_int); // F_int: [nsd*nnd]x[1]
}


void MFGP_Bal_EqT::Form_B_List(void)
{
		B1.Dimension(n_str, n_sd_x_n_en_displ);
		B3.Dimension(n_str, n_sd_x_n_en_displ);
		int dum=1;
		psi_lam.Dimension(dum, n_en_plast);
		B4.Dimension(dum, n_en_plast);
		
		/*
		fData_Pro_Displ->Set_B1(B1);
		fData_Pro_Displ->Set_B3(B3);
		fData_Pro_Plast->Set_psi_lam(psi_lam);
 		fData_Pro_Plast->Set_B4(B4); 
 		*/
 		Data_Pro_Displ.Set_B1(B1);
		Data_Pro_Displ.Set_B3(B3);
		Data_Pro_Plast.Set_psi_lam(psi_lam);
 		Data_Pro_Plast.Set_B4(B4); 
}


void MFGP_Bal_EqT::Form_C_List(GRAD_MRSSKStV *GRAD_MR_Plast)
{
		int i,j;
		Cuu1.Dimension(n_str,n_str); 
		Cuu2.Dimension(n_str,n_str);
		Culam1.Dimension(n_str,1); 
		Culam2.Dimension(n_str,1);
		Clamu1.Dimension(1,n_str); 
		Clamu2.Dimension(1,n_str);
		Clamlam1.Dimension(1); 
		Clamlam2.Dimension(1);
		
		/* retrive Cgep/fmoduli: 2D: [4]x[8]; 3D: [7]x[14] 
		*  and form the 8 C matrices
		*  pass the fModuli from the constitutive model */
		for (i=0; i<(n_str+1); ++i)
		{
		 for (j=0; j<(2*n_str+2); ++j)
		 {
		  if (i<n_str & j<n_str)
		  {
		   Cuu1(i,j)=moduli(i,j);
		  }
		  if (i<n_str & j>(n_str-1) & j<2*n_str)
		  {
		   Cuu2(i,j-n_str)=moduli(i,j);
		  }
		  if (i<n_str & j==2*n_str)
		  {
		   Culam1(i,0)=moduli(i,j);
		  }
		  if (i<n_str & j==(2*n_str+1))
		  {
		   Culam2(i,0)=moduli(i,j);
		  }
		  if (i==n_str & j<n_str)
		  {
		   Clamu1(0,j)=moduli(i,j);
		  }
		  if (i==n_str & j>(n_str-1) & j<2*n_str)
		  {
		   Clamu2(0,j-n_str)=moduli(i,j);
		  }
		  if (i==n_str & j==2*n_str)
		  {
		   Clamlam1(0,0)=moduli(i,j);
		  }
		  if (i==n_str & j==(2*n_str+1))
		  {
		   Clamlam2(0,0)=moduli(i,j);
		  }
		 } //inner for loop
		} //outer for loop
}



