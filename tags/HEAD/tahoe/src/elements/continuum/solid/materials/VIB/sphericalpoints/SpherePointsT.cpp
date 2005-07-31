/* $Id: SpherePointsT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (10/31/1997)                                          */
/* Base class for spherical point generators.                             */

#include "SpherePointsT.h"
#include <math.h>

const double Pi = acos(-1.0);

/*
* Constructor
*/
SpherePointsT::SpherePointsT(void): fQ(3)
{

}

/*
* Destructor
*/
SpherePointsT::~SpherePointsT() { }

/*
* List of jacobian determinants
*/
const dArrayT& SpherePointsT::Jacobians(void) const { return(fJacobians); }

/***********************************************************************
* Protected
***********************************************************************/

/*
* Orient points with
*/
void SpherePointsT::TransformPoints(double phi, double theta)
{
	/* convert to radians */
	phi   *= Pi/180.0;
	theta *= Pi/180.0;

	/* set tranformation tensor */
	fQ(0,0) = cos(phi);	
	fQ(1,0) = sin(phi);
	fQ(2,0) = 0.0;
	fQ(0,1) =-cos(theta)*sin(phi);
	fQ(1,1) = cos(theta)*cos(phi);
	fQ(2,1) =-sin(theta);
	fQ(0,2) =-sin(theta)*sin(phi);
	fQ(1,2) = cos(phi)*sin(theta);
	fQ(2,2) = cos(theta);

	dArrayT x, x_old(3);
	for (int i = 0; i < fPoints.MajorDim(); i++)
	{
		fPoints.RowAlias(i,x);
		x_old = x;
		fQ.Multx(x_old, x);
	}
}
