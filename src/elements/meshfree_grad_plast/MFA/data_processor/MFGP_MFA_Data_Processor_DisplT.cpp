// $Id: MFGP_MFA_Data_Processor_DisplT.cpp
#include "MFGP_MFA_Data_Processor_DisplT.h" 

using namespace Tahoe;

//---------------------------------------------------------------------

MFGP_MFA_Data_Processor_DisplT::MFGP_MFA_Data_Processor_DisplT() { };

//---------------------------------------------------------------------

MFGP_MFA_Data_Processor_DisplT::~MFGP_MFA_Data_Processor_DisplT() { };

MFGP_MFA_Data_Processor_DisplT::MFGP_MFA_Data_Processor_DisplT( dArray2DT &fdNdx, dArray2DT &fd3Ndx3 ) 
{
	Construct ( fdNdx, fd3Ndx3 );
}

//---------------------------------------------------------------------

void MFGP_MFA_Data_Processor_DisplT::Construct ( dArray2DT &fdNdx, dArray2DT &fd3Ndx3 )  
{
	dN = fdNdx;
	d3N = fd3Ndx3;
}


//---------------------------------------------------------------------
/* First Derivative of the Displacement Shape Function: [nsd] x [nnd] */ 
//fDphi has three components  
void MFGP_MFA_Data_Processor_DisplT::Set_B1( dMatrixT& B1 )
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


/* Laplacian of the Displacement Shape Function: [nstr] x [nnd] */ 
//fDDDphi has ten components; 
//void MFGP_MFA_Data_Processor_DisplT::Set_B3(MLSSolverGPT::SetShapeFunctions
//                        (const dArrayT& volume), dMatrixT& B3)
void MFGP_MFA_Data_Processor_DisplT::Set_B3( dMatrixT& B3 )
{
#if __option(extended_errorcheck)
	if (B3.Rows() != dSymMatrixT::NumValues(d3N.MajorDim()) ||
	    B3.Cols() != d3N.Length())
	    throw ExceptionT::kSizeMismatch;
#endif

	int nnd = d3N.MinorDim();
	double* pB3 = B3.Pointer();

	/* 1D */
	if (d3N.MajorDim() == 1)
	{
		const double* pNax = d3N(0);
		for (int i = 0; i < nnd; i++)
			*pB3++ = *pNax++;
	}
	/* 2D */
	else if (d3N.MajorDim() == 2)
	{
		const double* pNax = d3N(0);
		const double* pNay = d3N(1);
		const double* pNaxxy = d3N(5);
		const double* pNayyx = d3N(7);
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
		const double* pNax = d3N(0);
		const double* pNay = d3N(1);
		const double* pNaz = d3N(2);
		const double* pNaxxy = d3N(0); //components??
		const double* pNaxxz = d3N(1); //double check!!
		const double* pNayyx = d3N(2)
		const double* pNayyz = d3N(0);
		const double* pNazzx = d3N(1);
		const double* pNazzy = d3N(2);
		
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