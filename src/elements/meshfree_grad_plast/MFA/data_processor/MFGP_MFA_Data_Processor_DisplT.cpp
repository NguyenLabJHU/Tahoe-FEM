// $Id: MFGP_MFA_Data_Processor_DisplT.cpp
#include "MFGP_MFA_Data_Processor_DisplT.h" 

using namespace Tahoe;

/* constructor */
MFGP_MFA_Data_Processor_DisplT::MFGP_MFA_Data_Processor_DisplT() { };


/* destructor */
//MFGP_MFA_Data_Processor_DisplT::~MFGP_MFA_Data_Processor_DisplT() { };


/* initialize local variables */
void MFGP_MFA_Data_Processor_DisplT::Initialize(const dArray2DT& fdNdx, const dArray2DT& fd3Ndx3 )  
{
	dN = fdNdx;
	d3N = fd3Ndx3;
}


/* first derivative of the displacement shape function: [nsd] x [nnd] */  
void MFGP_MFA_Data_Processor_DisplT::Set_B1(dMatrixT& B1 )
{
#if __option(extended_errorcheck)
	if (B1.Rows() != dSymMatrixT::NumValues(dN.MajorDim()) ||
	    B1.Cols() != dN.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = dN.MinorDim();
	double* pB1 = B1.Pointer();

	/* 1D */
	if (dN.MajorDim() == 1)
	{
		const double* pNax = dN(0);
		for (int i = 0; i < nnd; i++)
			*pB1++ = *pNax++;
	}
	/* 2D */
	else if (dN.MajorDim() == 2)
	{
		const double* pNax = dN(0);
		const double* pNay = dN(1);
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
		const double* pNax = dN(0);
		const double* pNay = dN(1);
		const double* pNaz = dN(2);
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
void MFGP_MFA_Data_Processor_DisplT::Set_B3(dMatrixT& B3 )
{
#if __option(extended_errorcheck)
	if (B3.Rows() != dSymMatrixT::NumValues(sqrt(d3N.MajorDim())) ||
	    B3.Cols() != sqrt(d3N.Length()))
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = d3N.MinorDim();
	int nsd = sqrt(d3N.MajorDim());
	double* pB3 = B3.Pointer();

	/* 1D */
	if (nsd == 1)
	{
		const double* pNaxxx = d3N(0);
		for (int i = 0; i < nnd; i++)
			*pB3++ = *pNaxxx++;
	}
	/* 2D */
	else if (nsd == 2)
	{
		const double* pNaxxx = d3N(0);
		const double* pNayyx = d3N(1);
		const double* pNaxxy = d3N(2);
		const double* pNayyy = d3N(3);
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
		const double* pNaxxx = d3N(0);
		const double* pNayyx = d3N(1);
		const double* pNazzx = d3N(2);
		const double* pNaxxy = d3N(3); 
		const double* pNayyy = d3N(4); 
		const double* pNazzy = d3N(5);
		const double* pNaxxz = d3N(6);
		const double* pNayyz = d3N(7);
		const double* pNazzz = d3N(8);
		
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