/* $Id: D2OrthoMLS2DT.cpp,v 1.3 2002-10-20 22:49:42 paklein Exp $ */
/* created: paklein (10/17/1999)                                          */

#include "D2OrthoMLS2DT.h"

#include "ExceptionT.h"
#include "dSymMatrixT.h"

/* constructor */

using namespace Tahoe;

D2OrthoMLS2DT::D2OrthoMLS2DT(int complete):
	D2OrthoMLSSolverT(2, complete)
{
	/* supported bases - linear or quadratic */
	if (fComplete < 1 || fComplete > 2)
	{
		cout << "\n D2OrthoMLS2DT::D2OrthoMLS2DT: completeness of out of range {1,2}: ";
		cout << fComplete << endl;
		throw ExceptionT::kBadInputValue;
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* return the number of monomial terms for the given completeness */
int D2OrthoMLS2DT::NumberOfMonomials(int completeness) const
{
	switch (completeness)
	{
		case 0:			
			return 1;
		case 1:
			return 3;
		case 2:
			return 6;
		default:
			throw ExceptionT::kOutOfRange;
	}
	
	return 0;
}

/* return monomials evaluated at x */
void D2OrthoMLS2DT::SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (coords.Length() != fNumSD) throw ExceptionT::kGeneralFail;
	if (   p.Length() != NumberOfMonomials(fComplete)) throw ExceptionT::kSizeMismatch;
	if (Dp.MajorDim() != fNumSD ||
	    Dp.MinorDim() != p.Length()) throw ExceptionT::kSizeMismatch;
#endif

//NOTE: could do this for general completeness using
//      Outer with {f} and {f'} in each dimension

	double* pp = p.Pointer();
	double* px = Dp(0);
	double* py = Dp(1);
	
	double  x = coords[0];
	double  y = coords[1];

	switch (fComplete)
	{
		case 0: /* constant basis */

			p[0]  = 1.0;		

			px[0] = 0.0;
			py[0] = 0.0;
			break;	
		
		case 1: /* linear basis */
			
			*pp++ = 1.0;
			*pp++ = x;
			*pp   = y;

			*px++ = 0.0;
			*px++ = 1.0;
			*px   = 0.0;

			*py++ = 0.0;
			*py++ = 0.0;
			*py   = 1.0;			
			break;
			
		case 2: /* quadratic basis */
		
			*pp++ = 1.0;
			*pp++ = x;
			*pp++ = y;
			*pp++ = x*x;
			*pp++ = x*y;
			*pp   = y*y;

			*px++ = 0.0;
			*px++ = 1.0;
			*px++ = 0.0;
			*px++ = 2.0*x;
			*px++ = y;
			*px   = 0.0;

			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = 1.0;
			*py++ = 0.0;
			*py++ = x;
			*py   = 2.0*y;
			break;	

		default:
		
			throw ExceptionT::kOutOfRange;
	}
}

void D2OrthoMLS2DT::_SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp,
		dArray2DT& DDp)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (DDp.MajorDim() != dSymMatrixT::NumValues(fNumSD) ||
	    DDp.MinorDim() != p.Length()) throw ExceptionT::kSizeMismatch;
#endif

	/* set lower order derivatives */
	SetMonomials(coords, p, Dp);

	double* pxx = DDp(0);
	double* pyy = DDp(1);
	double* pxy = DDp(2);
	
	double  x = coords[0];
	double  y = coords[1];

	switch (fComplete)
	{
		case 0: /* constant basis */
		case 1: /* linear basis */
			
			DDp = 0.0;			
			break;
			
		case 2: /* quadratic basis */
		
			*pxx++ = 0.0;
			*pxx++ = 0.0;
			*pxx++ = 0.0;
			*pxx++ = 2.0;
			*pxx++ = 0.0;
			*pxx   = 0.0;

			*pyy++ = 0.0;
			*pyy++ = 0.0;
			*pyy++ = 0.0;
			*pyy++ = 0.0;
			*pyy++ = 0.0;
			*pyy   = 2.0;

			*pxy++ = 0.0;
			*pxy++ = 0.0;
			*pxy++ = 0.0;
			*pxy++ = 0.0;
			*pxy++ = 1.0;
			*pxy   = 0.0;

			break;	

		default:
		
			throw ExceptionT::kOutOfRange;
	}
}
