// $Id: MFGP_MFA_Data_Processor_DisplT.cpp
#include "MFGP_MFA_Data_Processor_DisplT.h" 

using namespace Tahoe;

//---------------------------------------------------------------------

MFGP_MFA_Data_Processor_DisplT::MFGP_MFA_Data_Processor_DisplT() { };

//---------------------------------------------------------------------

MFGP_MFA_Data_Processor_DisplT::~MFGP_MFA_Data_Processor_DisplT() { };


//---------------------------------------------------------------------
/* First Derivative of the Displacement Shape Function: [nsd] x [nnd] */ 
//fDphi has three components  
void MFGP_MFA_Data_Processor_DisplT::Set_B1(MLSSolverGPT::SetShapeFunctions
                        (const dArrayT& volume), dMatrixT& B1)
{
#if __option(extended_errorcheck)
	if (B1.Rows() != dSymMatrixT::NumValues(Dphi.MajorDim()) ||
	    B1.Cols() != Dphi.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = Dphi.MinorDim();
	double* pB1 = B1.Pointer();

	/* 1D */
	if (Dphi.MajorDim() == 1)
	{
		const double* pNax = Dphi(0);
		for (int i = 0; i < nnd; i++)
			*pB1++ = *pNax++;
	}
	/* 2D */
	else if (Dphi.MajorDim() == 2)
	{
		const double* pNax = Dphi(0);
		const double* pNay = Dphi(1);
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
		const double* pNax = Dphi(0);
		const double* pNay = Dphi(1);
		const double* pNaz = Dphi(2);
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


/* Laplacian of the Displacement Shape Function: [nstr] x [nnd] */ 
//fDDDphi has ten components; 
void MFGP_MFA_Data_Processor_DisplT::Set_B3(MLSSolverGPT::SetShapeFunctions
                        (const dArrayT& volume), dMatrixT& B3)
{
#if __option(extended_errorcheck)
	if (B3.Rows() != dSymMatrixT::NumValues(DDDphi.MajorDim()) ||
	    B3.Cols() != DDDphi.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DDDphi.MinorDim();
	double* pB3 = B3.Pointer();

	/* 1D */
	if (DDDphi.MajorDim() == 1)
	{
		const double* pNax = DDDphi(0);
		for (int i = 0; i < nnd; i++)
			*pB3++ = *pNax++;
	}
	/* 2D */
	else if (DDDphi.MajorDim() == 2)
	{
		const double* pNax = DDDphi(0);
		const double* pNay = DDDphi(1);
		const double* pNaxxy = DDDphi(5);
		const double* pNayyx = DDDphi(7);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB3++ = *pNax + (*pNaxxy);
			*pB3++ = 0.0;
			*pB3++ = *pNay + (*pNayyx);

			*pB3++ = 0.0;
			*pB3++ = *pNay + (*pNayyx);
			*pB3++ = *pNax + (*pNaxxy);
		}
	}
	/* 3D */
	else		
	{
		const double* pNax = DDDphi(0);
		const double* pNay = DDDphi(1);
		const double* pNaz = DDDphi(2);
		const double* pNaxxy = DDDphi(0); //components??
		const double* pNaxxz = DDDphi(1); //double check!!
		const double* pNayyx = DDDphi(2)
		const double* pNayyz = DDDphi(0);
		const double* pNazzx = DDDphi(1);
		const double* pNazzy = DDDphi(2)
		
		for (int i = 0; i < nnd; i++)
		{
			
			*pB3++ = *pNax + (*pNayyx) + (*pNazzx);
			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = *pNaz + (*pNaxxz) + (*pNayyz);
			*pB3++ = *pNay + (*pNaxxy) + (*pNazzy);

			*pB3++ = 0.0;
			*pB3++ = *pNay + (*pNaxxy) + (*pNazzy);
			*pB3++ = 0.0;
			*pB3++ = *pNaz + (*pNaxxz) + (*pNayyz);
			*pB3++ = 0.0;
			*pB3++ = *pNax + (*pNayyx) + (*pNazzx);

			*pB3++ = 0.0;
			*pB3++ = 0.0;
			*pB3++ = *pNaz++ + (*pNaxxz++) + (*pNayyz++);
			*pB3++ = *pNay++ + (*pNaxxy++) + (*pNazzy++);
			*pB3++ = *pNax++ + (*pNayyx++) + (*pNazzx++);
			*pB3++ = 0.0;
		}
	}
}