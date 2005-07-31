/* $Id: OrthoMLS3DT.cpp,v 1.1.1.1 2001-01-29 08:20:31 paklein Exp $ */
/* created: paklein (07/03/1998)                                          */

#include "OrthoMLS3DT.h"

/* constructor */
OrthoMLS3DT::OrthoMLS3DT(int complete): OrthoMLSSolverT(3, complete)
{
	/* supported bases - linear or quadratic */
	if (fComplete < 1 || fComplete > 3)
	{
		cout << "\n OrthoMLS3DT::OrthoMLS3DT: completeness of out of range {1,3}: ";
		cout << fComplete << endl;
		throw eBadInputValue;
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* return the number of monomial terms for the given completeness */
int OrthoMLS3DT::NumberOfMonomials(int completeness) const
{
	switch (completeness)
	{
		case 0:			
			return 1;
		case 1:
			return 4;
		case 2:
			return 10;
		case 3:
			return 20;
		default:
			throw eOutOfRange;
	}
	
	return 0;
}

/* return monomials evaluated at x */
void OrthoMLS3DT::SetMonomials(const dArrayT& coords, dArrayT& p, dArray2DT& Dp)
{
#if __option(extended_errorcheck)
	/* dimension checking */
	if (coords.Length() != fNumSD) throw eGeneralFail;
	if (   p.Length() != pow(fComplete+1,2)) throw eSizeMismatch;
	if (Dp.MajorDim() != fNumSD ||
	    Dp.MinorDim() != p.Length()) throw eSizeMismatch;
#endif

//NOTE: could do this for general completeness using
//      Outer with {f} and {f'} in each dimension

	double* pp = p.Pointer();
	double* px = Dp(0);
	double* py = Dp(1);
	double* pz = Dp(2);
	
	double  x = coords[0];
	double  y = coords[1];
	double  z = coords[2];
	
	switch (fComplete)
	{
		case 0: /* constant basis */

			p[0]  = 1.0;		

			px[0] = 0.0;
			py[0] = 0.0;
			pz[0] = 0.0;
			break;	
		
		case 1: /* linear basis */

			*pp++ = 1.0;
			*pp++ = x;
			*pp++ = y;
			*pp   = z;

			*px++ = 0.0;
			*px++ = 1.0;
			*px++ = 0.0;
			*px   = 0.0;

			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = 1.0;
			*py   = 0.0;

			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz   = 1.0;
			break;
			
		case 2: /* quadratic basis */

			*pp++ = 1.0;
			*pp++ = x;
			*pp++ = y;
			*pp++ = z;
			*pp++ = x*x;
			*pp++ = x*y;
			*pp++ = x*z;
			*pp++ = y*y;
			*pp++ = y*z;
			*pp   = z*z;

			*px++ = 0.0;
			*px++ = 1.0;
			*px++ = 0.0;
			*px++ = 0.0;
			*px++ = 2.0*x;
			*px++ = y;
			*px++ = z;
			*px++ = 0.0;
			*px++ = 0.0;
			*px   = 0.0;

			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = 1.0;
			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = x;
			*py++ = 0.0;
			*py++ = 2.0*y;
			*py++ = z;
			*py   = 0.0;

			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = 1.0;
			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = x;
			*pz++ = 0.0;
			*pz++ = y;
			*pz   = 2.0*z;
			break;	

		case 3: /* cubic basis */
		{
			double xx = x*x;
			double yy = y*y;
			double zz = z*z;

			double xy = x*y;
			double yz = z*y;
			double xz = x*z;

			*pp++ = 1.0;
			*pp++ = x;
			*pp++ = y;
			*pp++ = z;
			*pp++ = xx;
			*pp++ = xy;
			*pp++ = xz;
			*pp++ = yy;
			*pp++ = yz;
			*pp++ = zz;
			*pp++ = xx*x;
			*pp++ = xx*y;
			*pp++ = xx*z;
			*pp++ = x*yy;
			*pp++ = xy*z;
			*pp++ = x*zz;
			*pp++ = yy*y;
			*pp++ = yy*z;
			*pp++ = y*zz;
			*pp   = zz*z;

			*px++ = 0.0;
			*px++ = 1.0;
			*px++ = 0.0;
			*px++ = 0.0;
			*px++ = 2.0*x;
			*px++ = y;
			*px++ = z;
			*px++ = 0.0;
			*px++ = 0.0;
			*px++ = 0.0;
			*px++ = 3.0*xx;
			*px++ = 2.0*xy;
			*px++ = 2.0*xz;
			*px++ = yy;
			*px++ = yz;
			*px++ = zz;
			*px++ = 0.0;
			*px++ = 0.0;
			*px++ = 0.0;
			*px   = 0.0;

			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = 1.0;
			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = x;
			*py++ = 0.0;
			*py++ = 2.0*y;
			*py++ = z;
			*py++ = 0.0;
			*py++ = 0.0;
			*py++ = xx;
			*py++ = 0.0;
			*py++ = 2.0*xy;
			*py++ = xz;
			*py++ = 0.0;
			*py++ = 3.0*yy;
			*py++ = 2.0*yz;
			*py++ = zz;
			*py   = 0.0;

			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = 1.0;
			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = x;
			*pz++ = 0.0;
			*pz++ = y;
			*pz++ = 2.0*z;
			*pz++ = 0.0;
			*pz++ = 0.0;
			*pz++ = xx;
			*pz++ = 0.0;
			*pz++ = xy;
			*pz++ = 2.0*xz;
			*pz++ = 0.0;
			*pz++ = yy;
			*pz++ = 2.0*yz;
			*pz   = 3.0*zz;
			break;	
		}
		default:
		
			throw eOutOfRange;
	}
}
