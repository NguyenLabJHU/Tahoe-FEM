/* $Id: PolyBasis1DT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (12/11/1999)                                          */
/* base class for basis functions                                         */

#include "PolyBasis1DT.h"

/* constructor */
PolyBasis1DT::PolyBasis1DT(int complete):
	BasisT(complete, 1)
{
	if (fComplete < 0 || fComplete > 1)
	{
		cout << "\n PolyBasis1DT::PolyBasis1DT: completeness must be [0,1]" << endl;
		throw eBadInputValue;	
	}
}
	
/* return the number of basis functions */
int PolyBasis1DT::BasisDimension(void) const
{
	return fComplete + 1;
}

/* evaluate basis functions at coords */
void PolyBasis1DT::SetBasis(const dArray2DT& coords, int order)
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
				if (order > 1)
					fDDP[0] = 0.0;
			}
			break;
		}
		case 1: // linear basis
		{
			double* px   = coords.Pointer();
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
			break;
		}
	}
}
