/* $Id: LocalParabolaT.cpp,v 1.3 2002-10-20 22:39:04 paklein Exp $ */
/* created: paklein (01/28/1997)                                          */
/* LocalParabolaT.cpp                                                     */
/* For a successively fitted parabolic                                    */

#include "LocalParabolaT.h"
#include <math.h>
#include "ExceptionT.h"

/*
* Constructor
*/

using namespace Tahoe;

LocalParabolaT::LocalParabolaT(void)
{
	Reset();
}

/*
* Reset - clear all data
*/
void LocalParabolaT::Reset(void)
{
	/* status flags */
	fPointCount    = 0;
	fCurrentCoeffs = 0;

	/* data points */
	for (int i = 0; i < 3; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
	}
	
	/* coefficients y = a0 + a1 x + a2 x^2 */
	fa0 = 0.0;
	fa1 = 0.0;
	fa2 = 0.0;  		    	   	    		

	/* critical point values */
	fxmin = 0.0;  		    	   	    	
	fymin = 0.0;
}

/*
* Add point to fitting data, discarding the point in the current
* data farthest in y from the new point if at least 3 points
* have already been added.
*/
void LocalParabolaT::AddPoint(double x, double y)
{
	if (fPointCount < 3)
	{
		fx[fPointCount]   = x;
		fy[fPointCount]   = y;
	}
	else /* replace farthest y */
	{
		double dy0 = fabs(y - fy[0]);
		double dy1 = fabs(y - fy[1]);
		double dy2 = fabs(y - fy[2]);
		
		int dexmax = (dy0 > dy1) ?
					 ( (dy0 > dy2) ? 0 : 2 ) :
				     ( (dy1 > dy2) ? 1 : 2 );
									
		fx[dexmax] = x;
		fy[dexmax] = y;
	}

	/* status flags */
	fPointCount++;
	fCurrentCoeffs = 0;
}

/*
* Minima - return the critical point information through value.
* The functions return 0 if they are called before at least 3
* data points have been added, leaving value unmodified.
*/
int LocalParabolaT::CriticalPointValue(double& value) const
{
	if (fPointCount < 3)
		return 0;
	else
	{
		if (!fCurrentCoeffs)
			FitData();
	
		/* set return value */
		value = fymin;
	
		return 1;
	}
}

int LocalParabolaT::CriticalPointLocation(double& value) const
{
	if (fPointCount < 3)
		return 0;
	else
	{
		if (!fCurrentCoeffs)
			FitData();

		/* set return value */
		value = fxmin;
	
		return 1;
	}
}

int LocalParabolaT::CriticalPointConcavity(double& value) const
{
	if (fPointCount < 3)
		return 0;
	else
	{
		if (!fCurrentCoeffs)
			FitData();

		/* set return value */
		value = 2.0*fa2;
	
		return 1;
	}
}
	    	   	    	
/*
* Returning values
*/
double LocalParabolaT::operator()(double x) const
{
	if (fPointCount < 3)
		return 0.0;
	else
	{
		if (!fCurrentCoeffs)
			FitData();
	
		return fa0 + fa1*x + fa2*x*x;
	}
}
	
/**************************************************************************
* Private
**************************************************************************/

/*
* Compute the parabolic fit
*/
void LocalParabolaT::FitData(void) const
{
	/* make non-const local this */
	LocalParabolaT* localthis = (LocalParabolaT*) this;
	localthis->FitDataWrapper();
}

void LocalParabolaT::FitDataWrapper(void)
{
	/* checks */
	if (fPointCount < 3) throw ExceptionT::kGeneralFail;

	/* status flags */
	fCurrentCoeffs = 1;

	/* compute coefficients */
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21;

	double x1 = fx[0];
	double x2 = fx[1];
	double x3 = fx[2];

	double y1 = fy[0];
	double y2 = fy[1];
	double y3 = fy[2];

	z1 = -x1;
	z2 = x1*x1;
	z3 = -x2;
	z4 = x2*x2;
	z5 = x3*x3;
	z6 = x3*y1;
	z7 = x1*y2;
	z8 = -x3*y2;
	z9 = x2*y3;
	z10 = y3*z1;
	z11 = y1*z3;
	z12 = x2 + z1;
	z13 = x3 + z1;
	z14 = -y2*z2;
	z15 = y2*z2;
	z15 = x3*z15;
	z2 = y3*z2;
	z16 = z2*z3;
	z3 = x3 + z3;
	z17 = y1*z4;
	z18 = -y3*z4;
	z19 = x1*y3*z4;
	z20 = -y1*z5;
	z21 = x2*y1*z5;
	z5 = y2*z5;
	z1 = z1*z5;
	z4 = -z4*z6;
	z6 = z10 + z11 + z6 + z7 + z8 + z9;
	z7 = 1.0/z12;
	z8 = 1.0/z13;
	z3 = 1.0/z3;
	z2 = z14 + z17 + z18 + z2 + z20 + z5;
	z1 = z1 + z15 + z16 + z19 + z21 + z4;
	z3 = z3*z7*z8;
	z4 = z3*z6;
	z2 = z2*z3;
	z1 = z1*z3;
	
	//z1 = List(z1,z2,z4);
	fa0 = z1;
	fa1 = z2;
	fa2 = z4;

	/* critical point values */
	fxmin = -0.5*fa1/fa2;
	fymin = (*this)(fxmin);
}
