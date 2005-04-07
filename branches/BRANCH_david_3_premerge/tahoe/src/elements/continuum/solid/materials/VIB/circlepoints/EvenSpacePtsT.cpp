/* $Id: EvenSpacePtsT.cpp,v 1.6 2004-07-15 08:28:03 paklein Exp $ */
/* created: paklein (11/02/1997) */
#include "EvenSpacePtsT.h"

#include <math.h>
#include "toolboxConstants.h"
#include "ExceptionT.h"


using namespace Tahoe;

const double Pi = acos(-1.0);

/* constructor */
EvenSpacePtsT::EvenSpacePtsT(int n): fNtheta(n)
{
	/* number of integration points */
	if (fNtheta < 1) ExceptionT::BadInputValue("EvenSpacePtsT::EvenSpacePtsT");
	
	fPoints.Dimension(fNtheta,2);
	fJacobians.Dimension(fNtheta);
	
	/* all same weight */
	fJacobians = (2.0*Pi/fNtheta);
}

/*
* Generate points with the given orientation angle theta.
*/
const dArray2DT& EvenSpacePtsT::CirclePoints(double theta)
{
	/* generate direction vectors */
	double dtheta = (fNtheta == 2) ? Pi/2.0 : (2.0*Pi)/fNtheta;
	double angle  = theta - dtheta;
	dArrayT xsi;
	
	for (int i = 0; i < fNtheta; i++)
	{
		/* fetch vector */
		fPoints.RowAlias(i,xsi);
	
		/* orientation */
		angle += dtheta;
	
		/* components */
		xsi[0] = cos(angle);
		xsi[1] = sin(angle);
	}

	return fPoints;
}
