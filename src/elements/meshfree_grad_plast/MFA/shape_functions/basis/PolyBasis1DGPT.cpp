/* $Id: PolyBasis1DGPT.cpp,v 1.1 2004-06-22 23:17:48 kyonten Exp $ */
/* created: paklein (12/11/1999)                                          */
/* base class for basis functions                                         */

#include "PolyBasis1DGPT.h"

/* constructor */

using namespace Tahoe;

PolyBasis1DGPT::PolyBasis1DGPT(int complete):
	BasisGPT(complete, 1)
{
	if (fComplete < 0 || fComplete > 1)
	{
		cout << "\n PolyBasis1DGPT::PolyBasis1DGPT: completeness must be [0,1]" << endl;
		throw ExceptionT::kBadInputValue;	
	}
}
	
/* return the number of basis functions */
int PolyBasis1DGPT::BasisDimension(void) const
{
	return fComplete + 1;
}

/* evaluate basis functions at coords */
void PolyBasis1DGPT::SetBasis(const dArray2DT& coords, int order)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (coords.MinorDim() != fNumSD) throw ExceptionT::kGeneralFail;
	if (order > 2) throw ExceptionT::kOutOfRange;
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
				if (order > 1)
				{
					fDDP[0] = 0.0;
					if (order > 2) // kyonten (third derivative)
					{
						fDDDP[0] = 0.0;
					}
				} 
			}
			break;
		}
		case 1: // linear basis
		{
			const double* px = coords.Pointer();
			double* pP0  = fP(0);
			double* pP1  = fP(1);
			double* pDP0 = (fDP[0])(0);
			double* pDP1 = (fDP[0])(1);
			for (int i = 0; i < nnd; i++)
			{
				*pP0++ = 1.0;
				*pP1++ = *px++;
				if (order > 0)
				{
					*pDP0++ = 0.0;
					*pDP1++ = 1.0;
				}
			}
			
			if (order > 1) fDDP[0] = 0.0;
			if (order > 2) fDDDP[0] = 0.0; // kyonten
			break;
		}
	}
}
