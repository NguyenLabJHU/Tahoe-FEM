// $Id: MFGP_MFA_Data_Processor_PlastT.cpp
#include "MFGP_MFA_Data_Processor_PlastT.h"  

using namespace Tahoe;

//---------------------------------------------------------------------

MFGP_MFA_Data_Processor_PlastT::MFGP_MFA_Data_Processor_PlastT() { };

//---------------------------------------------------------------------

MFGP_MFA_Data_Processor_PlastT::~MFGP_MFA_Data_Processor_PlastT() { };


//---------------------------------------------------------------------
/* shape function of plastic multiplier */
// dof of lambda = 1 ??
// could be directly passed from the MLSSolverGP class??

void MFGP_MFA_Data_Processor_PlastT::Set_N(MLSSolverGPT::SetShapeFunctions
                        (const dArrayT& volume), dArrayT& N) 
{
	int nnd = phi.MinorDim();
	double* pN = N.Pointer();
	for (int i = 0; i < nnd; i++)
	  	*pN++ = *phi++;   			
}

//---------------------------------------------------------------------
/* Laplacian of the shape function of plastic multiplier */

void MFGP_MFA_Data_Processor_PlastT::Set_B4(MLSSolverGPT::SetShapeFunctions
                        (const dArrayT& volume), const B4)  
{
#if __option(extended_errorcheck)
	if (B4.Rows() != dSymMatrixT::NumValues(DDphi.MajorDim()) ||
	    B4.Cols() != DDphi.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = DDphi.MinorDim();
	double* pB4 = B4.Pointer();

	/* 1D */
	if (DDphi.MajorDim() == 1)
	{
		const double* pNax = DDphi(0);
		for (int i = 0; i < nnd; i++)
			*pB4++ = *pNax++;
	}
	/* 2D */
	else if (DDphi.MajorDim() == 2)
	{
		const double* pNaxx = DDphi(0);
		const double* pNayy = DDphi(1);
		for (int i = 0; i < nnd; i++)
		{
			
			*pB4++ = *pNaxx + (*pNayy);
		}
	}
	/* 3D */
	else		
	{
		const double* pNaxx = DDphi(0);
		const double* pNayy = DDphi(1);
		const double* pNazz = DDphi(2);
		
		for (int i = 0; i < nnd; i++)
		{
			
			*pB4++ = *pNaxx + (*pNayy) + (*pNazz);
		}
	}
}

