/* $Id: OrthoMLS2DT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (07/03/1998)                                          */

#include "OrthoMLS2DT.h"

/* constructor */
OrthoMLS2DT::OrthoMLS2DT(int complete): OrthoMLSSolverT(2, complete)
{
	/* supported bases - linear or quadratic */
	if (fComplete < 1 || fComplete > 3)
	{
		cout << "\n OrthoMLS2DT::OrthoMLS2DT: completeness of out of range {1,3}: ";
		cout << fComplete << endl;
		throw eBadInputValue;
	}
}

/*************************************************************************
* Protected
*************************************************************************/

/* return the number of monomial terms for the given completeness */
int OrthoMLS2DT::NumberOfMonomials(int completeness) const
{
	switch (completeness)
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

/* return monomials evaluated at x */
void OrthoMLS2DT::SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp)
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

		case 3: /* cubic basis */
		{
			double xx = x*x;
			double yy = y*y;
			double xy = x*y;

			*pp++ = 1.0;
			*pp++ = x;
			*pp++ = y;
			*pp++ = xx;
			*pp++ = xy;
			*pp++ = yy;
			*pp++ = xx*x;
			*pp++ = xx*y;
			*pp++ = x*yy;
			*pp   = yy*y;

			*px++ = 0.0;
			*px++ = 1.0;
			*px++ = 0.0;
			*px++ = 2.0*x;
			*px++ = y;
			*px++ = 0.0;
			*px++ = 3.0*xx;
			*px++ = 2.0*xy;
			*px++ = yy;
			*px   = 0.0;

			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = 1.0;
			*py++ = 0.0;
			*py++ = x;
			*py++ = 2.0*y;
			*py++ = 0.0;
			*py++ = xx;
			*py++ = 2.0*xy;
			*py   = 3.0*yy;
			break;	
		}
		default:
		
			throw eOutOfRange;
	}
}
