/* $Id: D2OrthoMLS1DT.cpp,v 1.1.1.1 2001-01-29 08:20:33 paklein Exp $ */
/* created: paklein (10/17/1999)                                          */

#include "D2OrthoMLS1DT.h"

#include "ExceptionCodes.h"
#include "dSymMatrixT.h"

/* constructor */
D2OrthoMLS1DT::D2OrthoMLS1DT(int complete):
	D2OrthoMLSSolverT(1, complete)
{
	/* supported bases - linear or quadratic */
	if (fComplete < 1 || fComplete > 3)
	{
		cout << "\n D2OrthoMLS1DT::D2OrthoMLS1DT: completeness of out of range {1,3}: ";
		cout << fComplete << endl;
		throw eBadInputValue;
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* return the number of monomial terms for the given completeness */
int D2OrthoMLS1DT::NumberOfMonomials(int completeness) const
{
	return completeness + 1;
}

/* return monomials evaluated at x */
void D2OrthoMLS1DT::SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (coords.Length() != fNumSD) throw eGeneralFail;
	if (   p.Length() != NumberOfMonomials(fComplete)) throw eSizeMismatch;
	if (Dp.MajorDim() != fNumSD ||
	    Dp.MinorDim() != p.Length()) throw eSizeMismatch;
#endif

//NOTE: could do this for general completeness using
//      Outer with {f} and {f'} in each dimension

	double* pp = p.Pointer();
	double* px = Dp(0);
	double   x = coords[0];

	switch (fComplete)
	{
		case 0: /* constant basis */

			 p[0] = 1.0;		
			px[0] = 0.0;
			break;	
		
		case 1: /* linear basis */
			
			*pp++ = 1.0;
			*pp++ = x;

			*px++ = 0.0;
			*px++ = 1.0;
			break;
			
		case 2: /* quadratic basis */
		
			*pp++ = 1.0;
			*pp++ = x;
			*pp++ = x*x;

			*px++ = 0.0;
			*px++ = 1.0;
			*px++ = 2.0*x;
			break;	

		default:
		
			throw eOutOfRange;
	}
}

void D2OrthoMLS1DT::_SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp,
		dArray2DT& DDp)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (DDp.MajorDim() != dSymMatrixT::NumValues(fNumSD) ||
	    DDp.MinorDim() != p.Length()) throw eSizeMismatch;
#endif

	/* set lower order derivatives */
	SetMonomials(coords, p, Dp);

	double* pxx = DDp(0);	
	double  x = coords[0];

	switch (fComplete)
	{
		case 0: /* constant basis */
		case 1: /* linear basis */
			
			DDp = 0.0;			
			break;
			
		case 2: /* quadratic basis */
		
			*pxx++ = 0.0;
			*pxx++ = 0.0;
			*pxx++ = 2.0;
			break;	

		default:
		
			throw eOutOfRange;
	}
}
