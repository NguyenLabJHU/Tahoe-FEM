// $Id: MFGP_MFA_Data_Processor_PlastT.cpp
#include "MFGP_MFA_Data_Processor_PlastT.h"  

using namespace Tahoe;

/* constructor */
MFGP_MFA_Data_Processor_PlastT::MFGP_MFA_Data_Processor_PlastT() { };


/* destructor */
//MFGP_MFA_Data_Processor_PlastT::~MFGP_MFA_Data_Processor_PlastT() { };


/* initialize local variables */
void MFGP_MFA_Data_Processor_PlastT::Initialize(const double* fN, const dArray2DT& fd2Ndx2 )  
{
	N = fN;
	d2N = fd2Ndx2;
}


/* shape function of plastic multiplier: [nnd] */ 
void MFGP_MFA_Data_Processor_PlastT::Set_psi_lam(dMatrixT& psi_lam) 
{
	int nnd = d2N.MinorDim();
	double* pphi = psi_lam.Pointer();
	for (int i = 0; i < nnd; i++)
	  	*pphi++ = *N++;		
}

/* laplacian of the shape function of plastic multiplier: [nnd] */
void MFGP_MFA_Data_Processor_PlastT::Set_B4(dMatrixT& B4)  
{
#if __option(extended_errorcheck)
	if (B4.Rows() != d2N.MajorDim() - d2N.MajorDim() ||
	    B4.Cols() != d2N.MinorDim())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = d2N.MinorDim();
	double* pB4 = B4.Pointer();

	/* 1D */
	if (d2N.MajorDim() == 1)
	{
		const double* pNaxx = d2N(0);
		for (int i = 0; i < nnd; i++)
			*pB4++ = *pNaxx++;
	}
	/* 2D */
	else if (d2N.MajorDim() == 2)
	{
		const double* pNaxx = d2N(0);
		const double* pNayy = d2N(1);
		for (int i = 0; i < nnd; i++)
		{
			*pB4++ = *pNaxx + (*pNayy);
		}
	}
	/* 3D */
	else		
	{
		const double* pNaxx = d2N(0);
		const double* pNayy = d2N(1);
		const double* pNazz = d2N(2);
		
		for (int i = 0; i < nnd; i++)
		{
			*pB4++ = *pNaxx + (*pNayy) + (*pNazz);
		}
	}
}

