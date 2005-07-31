/* $Id: PolyBasis3DT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (04/19/2000)                                          */

#include "PolyBasis3DT.h"

/* constructor */
PolyBasis3DT::PolyBasis3DT(int complete):
	BasisT(complete, 3)
{
	if (fComplete < 0 || fComplete > 1)
	{
		cout << "\n PolyBasis3DT::PolyBasis3DT: completeness must be [0,1]" << endl;
		throw eBadInputValue;	
	}
}
	
/* return the number of basis functions */
int PolyBasis3DT::BasisDimension(void) const
{
	switch (fComplete)
	{
		case 0:			
			return 1;
		case 1:
			return 4;
		case 2:
			return 10;
		default:
			throw eOutOfRange;
	}
	return 0;
}

/* evaluate basis functions at coords */
void PolyBasis3DT::SetBasis(const dArray2DT& coords, int order)
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
				fDP[2] = 0.0;
				if (order > 1)
				{
					fDDP[0] = 0.0;
					fDDP[1] = 0.0;
					fDDP[2] = 0.0;
					fDDP[3] = 0.0;
					fDDP[4] = 0.0;
					fDDP[5] = 0.0;
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
			double*  pP3 = fP(3);

			double* pD0P0 = (fDP[0])(0);
			double* pD0P1 = (fDP[0])(1);
			double* pD0P2 = (fDP[0])(2);
			double* pD0P3 = (fDP[0])(3);

			double* pD1P0 = (fDP[1])(0);
			double* pD1P1 = (fDP[1])(1);
			double* pD1P2 = (fDP[1])(2);
			double* pD1P3 = (fDP[1])(3);

			double* pD2P0 = (fDP[2])(0);
			double* pD2P1 = (fDP[2])(1);
			double* pD2P2 = (fDP[2])(2);
			double* pD2P3 = (fDP[2])(3);
			for (int i = 0; i < nnd; i++)
			{
				*pP0++ = 1.0;
				*pP1++ = *px++;
				*pP2++ = *px++;
				*pP3++ = *px++;

				if (order > 0)
				{
					*pD0P0++ = 0.0;
					*pD0P1++ = 1.0;
					*pD0P2++ = 0.0;
					*pD0P3++ = 0.0;

					*pD1P0++ = 0.0;
					*pD1P1++ = 0.0;
					*pD1P2++ = 1.0;
					*pD1P3++ = 0.0;

					*pD2P0++ = 0.0;
					*pD2P1++ = 0.0;
					*pD2P2++ = 0.0;
					*pD2P3++ = 1.0;
				}
			}
			if (order > 1)
			{
				fDDP[0] = 0.0;
				fDDP[1] = 0.0;
				fDDP[2] = 0.0;
				fDDP[3] = 0.0;
				fDDP[4] = 0.0;
				fDDP[5] = 0.0;
			}
			break;
		}
	}
}
