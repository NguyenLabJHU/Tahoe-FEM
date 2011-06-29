/* $Id: CirclePointsT.cpp,v 1.3 2004-07-15 08:28:03 paklein Exp $ */
/* created: paklein (11/02/1997) */
#include "CirclePointsT.h"
#include <math.h>

using namespace Tahoe;

const double Pi = acos(-1.0);

/*
* Constructor
*/
CirclePointsT::CirclePointsT(void): fQ(2)
{

}

/*
* Destructor
*/
CirclePointsT::~CirclePointsT() { }

/*
* List of jacobian determinants
*/
const dArrayT& CirclePointsT::Jacobians(void) const { return fJacobians; }

/***********************************************************************
* Protected
***********************************************************************/

/*
* Orient points with
*/
void CirclePointsT::TransformPoints(double theta)
{
	/* convert to radians */
	theta *= Pi/180.0;

	/* set tranformation tensor */
	fQ(0,0) = cos(theta);	fQ(0,1) =-sin(theta);
	fQ(1,0) = sin(theta);	fQ(1,1) = cos(theta);

	dArrayT x, x_old(2);
	for (int i = 0; i < fPoints.MajorDim(); i++)
	{
		fPoints.RowAlias(i,x);
		x_old = x;
		fQ.Multx(x_old, x);
	}
}
