/* $Id: EvenSpacePtsT.cpp,v 1.5 2004-06-17 07:40:52 paklein Exp $ */
/* created: paklein (11/02/1997)                                          */

#include "EvenSpacePtsT.h"

#include <math.h>
#include <iostream.h>

#include "toolboxConstants.h"
#include "ExceptionT.h"
#include "ifstreamT.h"


using namespace Tahoe;

const double Pi = acos(-1.0);

/*
* Constructor
*/
EvenSpacePtsT::EvenSpacePtsT(ifstreamT& in)
{
	/* number of integration points */
	in >> fNtheta;
	if (fNtheta < 1) throw ExceptionT::kBadInputValue;
	
	fPoints.Dimension(fNtheta,2);
	fJacobians.Dimension(fNtheta);
	
	/* all same weight */
	fJacobians = (2.0*Pi/fNtheta);
}

/*
* Print parameters.
*/
void EvenSpacePtsT::Print(ostream& out) const
{
	/* number of integration points */
	out << " Number of sampling points . . . . . . . . . . . = " << fNtheta << '\n';
}

void EvenSpacePtsT::PrintName(ostream& out) const
{
	out << "    Even spaced points\n";
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
