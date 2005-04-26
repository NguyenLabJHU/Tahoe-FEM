// $Id: Coupled_EqsT.cpp,v 1.1 2005-04-26 22:22:29 kyonten Exp $
#include "Coupled_EqsT.h" 

using namespace Tahoe;

/* constructor */
Coupled_EqsT::Coupled_EqsT(void) { }


/* destructor */
//Coupled_EqsT::~Coupled_EqsT(void) { }


void Coupled_EqsT::Initialize(int& curr_ip, D3MeshFreeShapeFunctionT* Shapes_displ, D3MeshFreeShapeFunctionT* Shapes_plast, 
							MFGPMaterialT* curr_mat, int& fTime_Step, double fdelta_t) 
{
	/* collect integers for dimensioning */
	n_en_displ  = Shapes_displ->Derivatives_U(curr_ip).MinorDim(); //dN: [nsd]x[nnd]
	n_en_plast  = Shapes_plast->Derivatives_U(curr_ip).MinorDim();
	n_sd    	= Shapes_displ->Derivatives_U(curr_ip).MajorDim();
	n_str = dSymMatrixT::NumValues(n_sd);	
	n_sd_x_n_en_displ = n_sd * n_en_displ;
	n_sd_x_n_en_plast = n_sd * n_en_plast;

	delta_t = fdelta_t;
 	
	CurrMat = curr_mat;
	N    = Shapes_displ->IPShapeU(curr_ip);
	DN   = Shapes_displ->Derivatives_U(curr_ip);
	DDN  = Shapes_displ->DDerivatives_U(curr_ip);
	DDDN = Shapes_displ->DDDerivatives_U(curr_ip);
	
	Form_C_List();
	Form_B_List();
}

//---------------------------------------------------------------------

void Coupled_EqsT::Form_KUU_KULam(dMatrixT& Kuu, dMatrixT& Kulam)  
{		
	#if __option(extended_errorcheck)
		if (Kuu.Rows() != B1.Cols() || Kulam.Rows() != B1.Cols() ||
		    Kulam.Cols() != psi_lam.Cols())
			throw ExceptionT::kSizeMismatch;
	#endif
	
	dMatrixT Ktemp1(Kuu.Rows(), Kuu.Cols());
	Kuu.MultATBC(B1, Cuu1, B1);
	Ktemp1.MultATBC(B1, Cuu2, B3);
	Kuu += Ktemp1; // Kuu: [nsd*nnd]x[nsd*nnd]
	
	dMatrixT Ktemp2(Kulam.Rows(), Kulam.Cols());
	Kulam.MultATBC(B1, Culam1, psi_lam);
	Ktemp2.MultATBC(B1, Culam2, B4);
	Kulam += Ktemp2;	// Kulambda: [nsd*nnd]x[nnd] 
}

void Coupled_EqsT::Form_FU_int(dArrayT& Fu_int) 
{
	
	#if __option(extended_errorcheck)
		if (Fu_int.Length() != B1.Cols()) throw ExceptionT::kSizeMismatch;
	#endif
	
	dMatrixT B1_T(n_sd_x_n_en_displ, n_str);
	dArrayT stress_sym(n_str);
	B1_T.Transpose(B1);
	dSymMatrixT stress = CurrMat->s_ij();
	/* take out the symmetric part of the stress */
	for (int ij = 0; ij < n_str; ij++)
	{
		int i, j;
		dSymMatrixT::ExpandIndex(n_sd, ij, i, j);
		stress_sym[ij] = stress(i,j);
	} 
	
	B1_T.MultTx(stress_sym, Fu_int); // Fu_int: [nsd*nnd]x[1]
}


void Coupled_EqsT::Form_KLamU_KLamLam(dMatrixT& Klamu, dMatrixT& Klamlam)  
{
	#if __option(extended_errorcheck)
		if (Klamlam.Rows() != psi_lam.Cols() || Klamu.Rows() != psi_lam.Cols() ||
		    Klamu.Cols() != B1.Cols())
			throw ExceptionT::kSizeMismatch;
	#endif
	
	dMatrixT Ktemp1(Klamlam.Rows(), Klamlam.Cols());
	Klamlam.MultATBC(psi_lam, Clamlam1, psi_lam);
	Ktemp1.MultATBC(psi_lam, Clamlam2, B4);
	Klamlam += Ktemp1;	//Klamlam: [nnd]x[nnd]
	
	dMatrixT Ktemp2(Klamu.Rows(), Klamu.Cols());
	Klamu.MultATBC(psi_lam, Clamu1, B1);
	Ktemp2.MultATBC(psi_lam, Clamu2, B3);
	Klamu += Ktemp2;	// Klamu :[nnd]x[nsd*nnd]
}

//---------------------------------------------------------------------

void Coupled_EqsT::Form_FLambda_int(dArrayT& Flambda_int) 
{
	#if __option(extended_errorcheck)
		if (Flambda_int.Length() != B4.Cols()) throw ExceptionT::kSizeMismatch;
	#endif
		
	//pass column of phi_lam to F_int
	double yield = CurrMat->YieldF();
	Flambda_int = psi_lam[0];
	Flambda_int *= yield; // Flambda_int: [nnd]x[1]
	//updated failure function passed from the constitutive model
}

void Coupled_EqsT::Form_B_List(void)
{
	/* dimension */
	B1.Dimension(n_str, n_sd_x_n_en_displ);
	B3.Dimension(n_str, n_sd_x_n_en_displ);
	int dum=1;
	psi_lam.Dimension(dum, n_en_plast);
	B4.Dimension(dum, n_en_plast); 
}


void Coupled_EqsT::Form_C_List()
{
	/* dimensions */
	Cuu1.Dimension(n_str,n_str); 
	Cuu2.Dimension(n_str,n_str);
	Culam1.Dimension(n_str,1); 
	Culam2.Dimension(n_str,1);
	Clamu1.Dimension(1,n_str); 
	Clamu2.Dimension(1,n_str);
	Clamlam1.Dimension(1,1); 
	Clamlam2.Dimension(1,1);
		
	/* C matrices */
	Cuu1 = CurrMat->c_UU1_ijkl(); 
	Cuu2 = CurrMat->c_UU2_ijkl();
	Culam1 = CurrMat->c_ULam1_ij(); 
	Culam2 = CurrMat->c_ULam2_ij();
	Clamu1 = CurrMat->c_LamU1_ij(); 
	Clamu2 = CurrMat->c_LamU2_ij();
	Clamlam1 = CurrMat->c_LamLam1(); 
	Clamlam2 = CurrMat->c_LamLam2();
		
}

/* first derivative of the displacement shape function: [nsd] x [nnd] */  
void Coupled_EqsT::Set_B1(dMatrixT& B1 )
{
#if __option(extended_errorcheck)
	if (B1.Rows() != dSymMatrixT::NumValues(DN.MajorDim()) ||
	    B1.Cols() != DN.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DN.MinorDim();
	double* pB1 = B1.Pointer();

	/* 1D */
	if (DN.MajorDim() == 1)
	{
		const double* pNax = DN(0);
		for (int i = 0; i < nnd; i++)
			*pB1++ = *pNax++;
	}
	/* 2D */
	else if (DN.MajorDim() == 2)
	{
		const double* pNax = DN(0);
		const double* pNay = DN(1);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB1++ = *pNax;
			*pB1++ = 0.0;
			*pB1++ = *pNay;

			*pB1++ = 0.0;
			*pB1++ = *pNay++;
			*pB1++ = *pNax++;
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DN(0);
		const double* pNay = DN(1);
		const double* pNaz = DN(2);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB1++ = *pNax;
			*pB1++ = 0.0;
			*pB1++ = 0.0;
			*pB1++ = 0.0;
			*pB1++ = *pNaz;
			*pB1++ = *pNay;

			*pB1++ = 0.0;
			*pB1++ = *pNay;
			*pB1++ = 0.0;
			*pB1++ = *pNaz;
			*pB1++ = 0.0;
			*pB1++ = *pNax;

			*pB1++ = 0.0;
			*pB1++ = 0.0;
			*pB1++ = *pNaz++;
			*pB1++ = *pNay++;
			*pB1++ = *pNax++;
			*pB1++ = 0.0;
		}
	}
}


/* laplacian of the displacement shape function: [nsd*nsd] x [nnd] */  
void Coupled_EqsT::Set_B3(dMatrixT& B3)
{
#if __option(extended_errorcheck)
//	if (B3.Rows() != dSymMatrixT::NumValues(sqrt(d3N.MajorDim())) ||
//	    B3.Cols() != sqrt(d3N.MajorDim())*d3N.MinorDim())
//	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DDDN.MinorDim();
	int nsd = 3;
	if (DDDN.MajorDim() == 4)
		nsd = 2;
	else if (DDDN.MajorDim() == 1)
		nsd = 1;
	double* pB3 = B3.Pointer();

	/* 1D */
	if (nsd == 1)
	{
		const double* pNaxxx = DDDN(0);
		for (int i = 0; i < nnd; i++)
			*pB3++ = *pNaxxx++;
	}
	/* 2D */
	else if (nsd == 2)
	{
		const double* pNaxxx = DDDN(0);
		const double* pNayyx = DDDN(1);
		const double* pNaxxy = DDDN(2);
		const double* pNayyy = DDDN(3);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB3++ = *pNaxxx + (*pNayyx);
			*pB3++ = 0.0;
			*pB3++ = *pNayyy + (*pNaxxy);

			*pB3++ = 0.0;
			*pB3++ = *pNayyy + (*pNaxxy);
			*pB3++ = *pNaxxx + (*pNayyx);
		}
	}
	/* 3D */
	else		
	{
		const double* pNaxxx = DDDN(0);
		const double* pNayyx = DDDN(1);
		const double* pNazzx = DDDN(2);
		const double* pNaxxy = DDDN(3); 
		const double* pNayyy = DDDN(4); 
		const double* pNazzy = DDDN(5);
		const double* pNaxxz = DDDN(6);
		const double* pNayyz = DDDN(7);
		const double* pNazzz = DDDN(8);
		
		for (int i = 0; i < nnd; i++)
		{
			
			*pB3++ = *pNaxxx + (*pNayyx) + (*pNazzx);
			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = *pNazzz + (*pNaxxz) + (*pNayyz);
			*pB3++ = *pNayyy + (*pNaxxy) + (*pNazzy);

			*pB3++ = 0.0;
			*pB3++ = *pNayyy + (*pNaxxy) + (*pNazzy);
			*pB3++ = 0.0;
			*pB3++ = *pNazzz + (*pNaxxz) + (*pNayyz);
			*pB3++ = 0.0;
			*pB3++ = *pNaxxx + (*pNayyx) + (*pNazzx);

			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = *pNazzz++ + (*pNaxxz++) + (*pNayyz++);
			*pB3++ = *pNayyy++ + (*pNaxxy++) + (*pNazzy++);
			*pB3++ = *pNaxxx++ + (*pNayyx++) + (*pNazzx++);
			*pB3++ = 0.0;
		}
	}
}

/* shape function of plastic multiplier: [1]x[nnd] */ 
void Coupled_EqsT::Set_psi_lam(dMatrixT& psi_lam) 
{
	int nnd = DN.MinorDim();
	double* pphi = psi_lam.Pointer();
	for (int i = 0; i < nnd; i++)
	  	*pphi++ = *N++;
		
}

/* laplacian of the shape function of plastic multiplier: [1]x[nnd] */
void Coupled_EqsT::Set_B4(dMatrixT& B4)  
{
#if __option(extended_errorcheck)
	if (B4.Rows() != DDN.MajorDim() - DDN.MajorDim() ||
	    B4.Cols() != DDN.MinorDim())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DDN.MinorDim();
	double* pB4 = B4.Pointer();

	/* 1D */
	if (DDN.MajorDim() == 1)
	{
		const double* pNaxx = DDN(0);
		for (int i = 0; i < nnd; i++)
			*pB4++ = *pNaxx++;
	}
	/* 2D */
	else if (DDN.MajorDim() == 2)
	{
		const double* pNaxx = DDN(0);
		const double* pNayy = DDN(1);
		for (int i = 0; i < nnd; i++)
		{
			*pB4++ = *pNaxx + (*pNayy);
		}
	}
	/* 3D */
	else		
	{
		const double* pNaxx = DDN(0);
		const double* pNayy = DDN(1);
		const double* pNazz = DDN(2);
		
		for (int i = 0; i < nnd; i++)
		{
			*pB4++ = *pNaxx + (*pNayy) + (*pNazz);
		}
	}
}

