/* $Id: PolyBasis2DT.cpp,v 1.8 2005-04-22 00:33:40 paklein Exp $ */
/* created: paklein (12/13/1999) */
#include "PolyBasis2DT.h"

using namespace Tahoe;

/* constructor */
PolyBasis2DT::PolyBasis2DT(int complete, bool cross_terms):
	BasisT(complete, 2),
	fCrossTerms(cross_terms)
{
	/* check */
	if (fComplete < 0 || fComplete > 1)
		ExceptionT::OutOfRange("PolyBasis2DT::PolyBasis2DT",
			"completeness must be [0,1]: %d", complete);
}
	
/* return the number of basis functions */
int PolyBasis2DT::BasisDimension(void) const
{
	switch (fComplete)
	{
		case 0:			
			return 1;
		case 1:
			if (!fCrossTerms)
				return 3;
			else
				return 4;
		case 2:
			if (!fCrossTerms)
				return 6;
			else
				return 9;
		case 3:
			if (!fCrossTerms)
				return 10;
			else
				return 16;
		default:
			ExceptionT::OutOfRange("PolyBasis2DT::BasisDimension");
	}
	
	return 0;
}

/* evaluate basis functions at coords */
void PolyBasis2DT::SetBasis(const dArray2DT& coords, int order)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	const char caller[] = "PolyBasis2DT::SetBasis";
	if (coords.MinorDim() != fNumSD) ExceptionT::GeneralFail(caller);
	if (order > 3) ExceptionT::OutOfRange(caller); //kyonten (order increase to 3)
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
					if (order > 2) // kyonten
					{
						fDDDP[0] = 0.0;
						fDDDP[1] = 0.0;
						fDDDP[2] = 0.0;
						fDDDP[3] = 0.0;
					}
				}
			}
			break;
		}
		case 1: // linear basis
		{
			const double* px = coords.Pointer();
			double*  pP0 = fP(0);
			double*  pP1 = fP(1);
			double*  pP2 = fP(2);
			double*  pP3 = NULL;
			
			double* pD0P0 = (fDP[0])(0);
			double* pD0P1 = (fDP[0])(1);
			double* pD0P2 = (fDP[0])(2);
			double* pD0P3 = NULL;

			double* pD1P0 = (fDP[1])(0);
			double* pD1P1 = (fDP[1])(1);
			double* pD1P2 = (fDP[1])(2);
			double* pD1P3 = NULL;

			if (fCrossTerms) {
				pP3   = (fCrossTerms) ?       fP(3) : NULL;
				pD0P3 = (fCrossTerms) ? (fDP[0])(3) : NULL;
				pD1P3 = (fCrossTerms) ? (fDP[1])(3) : NULL;
			}

			for (int i = 0; i < nnd; i++)
			{
				double x = *px++;
				double y = *px++;
			
				*pP0++ = 1.0;
				*pP1++ = x;
				*pP2++ = y;
				if (fCrossTerms) *pP3++ = x*y;
				
				if (order > 0)
				{
					*pD0P0++ = 0.0;
					*pD0P1++ = 1.0;
					*pD0P2++ = 0.0;

					*pD1P0++ = 0.0;
					*pD1P1++ = 0.0;
					*pD1P2++ = 1.0;
					
					if (fCrossTerms) {
						*pD0P3++ = y;
						*pD1P3++ = x;
					}

					if (order > 1)
					{
						fDDP[0] = 0.0;
						fDDP[1] = 0.0;
						fDDP[2] = 0.0;
						
						if (fCrossTerms)
							(fDDP[2])(3,i) = 1.0;
						
						if (order > 2)
						{
							fDDDP[0] = 0.0;
							fDDDP[1] = 0.0;
							fDDDP[2] = 0.0;
							fDDDP[3] = 0.0;
						}
					}
				}
			}			
			break;
		}
	}
}
