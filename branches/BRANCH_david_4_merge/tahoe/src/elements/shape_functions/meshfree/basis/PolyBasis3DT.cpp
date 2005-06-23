/* $Id: PolyBasis3DT.cpp,v 1.8 2005-04-22 00:33:40 paklein Exp $ */
/* created: paklein (04/19/2000) */
#include "PolyBasis3DT.h"

using namespace Tahoe;

/* constructor */
PolyBasis3DT::PolyBasis3DT(int complete, bool cross_terms):
	BasisT(complete, 3),
	fCrossTerms(cross_terms)
{
	/* check */
	if (fComplete < 0 || fComplete > 1)
		ExceptionT::OutOfRange("PolyBasis3DT::PolyBasis3DT",
			"completeness must be [0,1]: %d", complete);
}
	
/* return the number of basis functions */
int PolyBasis3DT::BasisDimension(void) const
{
	const char caller[] = "PolyBasis3DT::BasisDimension";
	switch (fComplete)
	{
		case 0:			
			return 1;
		case 1:
			if (!fCrossTerms)
				return 4;
			else
				return 7;
		case 2:
			if (!fCrossTerms)
				return 10;
			else
				return 20;
		case 3:
			if (fCrossTerms) ExceptionT::GeneralFail(caller);
			return 10; 
		default:
			ExceptionT::OutOfRange(caller);
	}
	return 0;
}

/* evaluate basis functions at coords */
void PolyBasis3DT::SetBasis(const dArray2DT& coords, int order)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	const char caller[] = "PolyBasis3DT::SetBasis";
	if (coords.MinorDim() != fNumSD) ExceptionT::GeneralFail(caller);
	if (order > 3) ExceptionT::OutOfRange(caller); //kyonten (order increased to 3)
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
					if (order > 2) // kyonten
					{
						fDDDP[0] = 0.0;
						fDDDP[1] = 0.0;
						fDDDP[2] = 0.0;
						fDDDP[3] = 0.0;
						fDDDP[4] = 0.0;
						fDDDP[5] = 0.0;
						fDDDP[6] = 0.0;
						fDDDP[7] = 0.0;
						fDDDP[8] = 0.0;
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
			double*  pP3 = fP(3);
			double*  pP4 = NULL;
			double*  pP5 = NULL;
			double*  pP6 = NULL;

			double* pD0P0 = (fDP[0])(0);
			double* pD0P1 = (fDP[0])(1);
			double* pD0P2 = (fDP[0])(2);
			double* pD0P3 = (fDP[0])(3);
			double* pD0P4 = NULL;
			double* pD0P5 = NULL;
			double* pD0P6 = NULL;

			double* pD1P0 = (fDP[1])(0);
			double* pD1P1 = (fDP[1])(1);
			double* pD1P2 = (fDP[1])(2);
			double* pD1P3 = (fDP[1])(3);
			double* pD1P4 = NULL;
			double* pD1P5 = NULL;
			double* pD1P6 = NULL;

			double* pD2P0 = (fDP[2])(0);
			double* pD2P1 = (fDP[2])(1);
			double* pD2P2 = (fDP[2])(2);
			double* pD2P3 = (fDP[2])(3);
			double* pD2P4 = NULL;
			double* pD2P5 = NULL;
			double* pD2P6 = NULL;

			if (fCrossTerms) {
				pP4 = fP(4);
				pP5 = fP(5);
				pP6 = fP(6);

				pD0P4 = (fDP[0])(4);
				pD0P5 = (fDP[0])(5);
				pD0P6 = (fDP[0])(6);

				pD1P4 = (fDP[1])(4);
				pD1P5 = (fDP[1])(5);
				pD1P6 = (fDP[1])(6);

				pD2P4 = (fDP[2])(4);
				pD2P5 = (fDP[2])(5);
				pD2P6 = (fDP[2])(6);			
			}

			for (int i = 0; i < nnd; i++)
			{
				double x = *px++;
				double y = *px++;
				double z = *px++;
			
				*pP0++ = 1.0;
				*pP1++ = x;
				*pP2++ = y;
				*pP3++ = z;
				if (fCrossTerms) {
					*pP4++ = y*z;
					*pP5++ = x*z;
					*pP6++ = x*y;
				}

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
					
					if (fCrossTerms) {
						*pD0P4++ = 0.0;
						*pD0P5++ = z;
						*pD0P6++ = y;

						*pD1P4++ = z;
						*pD1P5++ = 0.0;
						*pD1P6++ = x;

						*pD2P4++ = y;
						*pD2P5++ = x;
						*pD2P6++ = 0.0;
					}
					
					if (order > 1)
					{
						fDDP[0] = 0.0;
						fDDP[1] = 0.0;
						fDDP[2] = 0.0;
						fDDP[3] = 0.0;
						fDDP[4] = 0.0;
						fDDP[5] = 0.0;

						if (fCrossTerms) {
							(fDDP[3])(4,i) = 1.0;
							(fDDP[4])(5,i) = 1.0;
							(fDDP[5])(6,i) = 1.0;						
						}

						if (order > 2) // kyonten
						{
							fDDDP[0] = 0.0;
							fDDDP[1] = 0.0;
							fDDDP[2] = 0.0;
							fDDDP[3] = 0.0;
							fDDDP[4] = 0.0;
							fDDDP[5] = 0.0;
							fDDDP[6] = 0.0;
							fDDDP[7] = 0.0;
							fDDDP[8] = 0.0;
						}
					}	
				}
			}
			break;
		}
	}
}
