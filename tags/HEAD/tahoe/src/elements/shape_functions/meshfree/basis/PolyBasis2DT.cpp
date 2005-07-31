/* $Id: PolyBasis2DT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (12/13/1999)                                          */

#include "PolyBasis2DT.h"

/* constructor */
PolyBasis2DT::PolyBasis2DT(int complete):
	BasisT(complete, 2)
{
	if (fComplete < 0 || fComplete > 1)
	{
		cout << "\n PolyBasis2DT::PolyBasis2DT: completeness must be [0,1]" << endl;
		throw eBadInputValue;	
	}
}
	
/* return the number of basis functions */
int PolyBasis2DT::BasisDimension(void) const
{
	switch (fComplete)
	{
		case 0:			
			return 1;
		case 1:
			return 3;
		case 2:
			return 6;
		case 3:
			return 10;
		default:
			throw eOutOfRange;
	}
	
	return 0;
}

/* evaluate basis functions at coords */
void PolyBasis2DT::SetBasis(const dArray2DT& coords, int order)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (coords.MinorDim() != fNumSD) throw eGeneralFail;
	if (order > 2) throw eOutOfRange;
#endif

	/* dimensions */
	int nnd = coords.MajorDim();

	/* dimension work space */
	Dimension(nnd);

	switch (fComplete)
	{
		case 0: // constant basis
		{
			fP = 1.0;
			if (order > 0)
			{
				fDP[0] = 0.0;
				fDP[1] = 0.0;
				if (order > 1)
				{
					fDDP[0] = 0.0;
					fDDP[1] = 0.0;
					fDDP[2] = 0.0;
				}
			}
			break;
		}
		case 1: // linear basis
		{
			double*   px = coords.Pointer();
			double*  pP0 = fP(0);
			double*  pP1 = fP(1);
			double*  pP2 = fP(2);
			
			double* pD0P0 = (fDP[0])(0);
			double* pD0P1 = (fDP[0])(1);
			double* pD0P2 = (fDP[0])(2);

			double* pD1P0 = (fDP[1])(0);
			double* pD1P1 = (fDP[1])(1);
			double* pD1P2 = (fDP[1])(2);
			for (int i = 0; i < nnd; i++)
			{
				*pP0++ = 1.0;
				*pP1++ = *px++;
				*pP2++ = *px++;
				
				if (order > 0)
				{
					*pD0P0++ = 0.0;
					*pD0P1++ = 1.0;
					*pD0P2++ = 0.0;

					*pD1P0++ = 0.0;
					*pD1P1++ = 0.0;
					*pD1P2++ = 1.0;
				}
			}
			
			if (order > 1)
			{
				fDDP[0] = 0.0;
				fDDP[1] = 0.0;
				fDDP[2] = 0.0;
			}
			break;
		}
	}
}
