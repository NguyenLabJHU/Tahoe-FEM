/* $Id: GaussPtsT.cpp,v 1.5 2004-07-15 08:28:03 paklein Exp $ */
/* created: paklein (11/02/1997) */
#include "GaussPtsT.h"
#include <math.h>
#include "toolboxConstants.h"
#include "ExceptionT.h"

using namespace Tahoe;

const double Pi = acos(-1.0);

/* constructor */
GaussPtsT::GaussPtsT(int n): fN(n)
{

}

/*
* Generate points with the given orientation angle theta.
*/
const dArray2DT& GaussPtsT::CirclePoints(double theta)
{	
	/* set jacobians */
	SetJacobians(fN);

	/* set coordinates */
	SetCoords(fN);
	
	/* reorient points */
	TransformPoints(theta);
	
	return fPoints;
}

/*
* Returns the correct data pointer for the specified number of
* integration points
*/
void GaussPtsT::SetCoords(int numint)
{
	/* parent domain points*/
	double p9[]=
	{-0.96816023950762608984,
	 -0.83603110732663579430,
	 -0.61337143270059039731,
	 -0.32425342340380892904,
	  0.0,
	  0.32425342340380892904,
	  0.61337143270059039731,
	  0.83603110732663579430,
	  0.96816023950762608984};

	double p10[]=
	{-0.97390652851717172008,
	 -0.86506336668898451073,
	 -0.67940956829902440623,
	 -0.43339539412924719080,
	 -0.14887433898163121089,
	  0.14887433898163121089,
	  0.43339539412924719080,
	  0.67940956829902440623,
	  0.86506336668898451073,
	  0.97390652851717172008};
	
	double *p;
	switch (numint)
	{
		case 9:

			p = p9;
			break;

		case 10:

			p = p10;
			break;

		default:
			ExceptionT::GeneralFail("GaussPtsT::SetCoords", "unrecognized Gauss rule %d", numint);
	}

	/* calculate directions */	
	fPoints.Dimension(numint,2);
	for (int i = 0; i < numint; i++)
	{
		double *xsi = fPoints(i);
		
		/* set direction cosines */
		xsi[0] = cos(Pi*p[i]);
		xsi[1] = sin(Pi*p[i]);	
	}
}

void GaussPtsT::SetJacobians(int numint)
{
	double p9[]=
{0.08127438836157441197,
0.18064816069485740406,
0.26061069640293546232,
0.31234707704000284007,
0.33023935500125976317,
0.31234707704000284007,
0.26061069640293546232,
0.18064816069485740406,
0.08127438836157441197};

	double p10[]=
{0.06667134430868813759,
0.14945134915058059315,
0.21908636251598204400,
0.26926671930999635509,
0.29552422471475287017,
0.29552422471475287017,
0.26926671930999635509,
0.21908636251598204400,
0.14945134915058059315,
0.06667134430868813759};

	double *p;
	switch (numint)
	{
		case 9:

			p = p9;
			break;

		case 10:

			p = p10;
			break;

		default:
			ExceptionT::GeneralFail("GaussPtsT::SetCoords", "unrecognized Gauss rule %d", numint);
	}
	
	/* temp vector */
	dArrayT temp(numint,p);
	
	/* copy in */
	fJacobians.Dimension(numint);
	fJacobians.SetToScaled(Pi,temp);
}
