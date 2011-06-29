/* $Id: EvenSpacePtsT.cpp,v 1.7 2006-08-18 18:45:11 tdnguye Exp $ */
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
	fAngles.Dimension(fNtheta);
	fJacobians.Dimension(fNtheta);
	
	/* all same weight */
	fJacobians = (2.0*Pi/fNtheta);
}

const dArrayT& EvenSpacePtsT::Jacobians(const double theta, const C1FunctionT* func) 
{
	/* generate direction vectors */
	double dtheta = (fNtheta == 2) ? Pi/2.0 : (2.0*Pi)/fNtheta;
	double angle  = theta - dtheta;
	
	for (int i = 0; i < fNtheta; i++)
	{
		/* orientation */
		angle += dtheta;
		double D = func->Function(angle);
		
		/* components */
		fJacobians[i] = D*(2.0*Pi/fNtheta);
	}

	return fJacobians;
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

const dArrayT& EvenSpacePtsT::CircleAngles(double theta)
{
	/* generate direction vectors */
	double dtheta = (fNtheta == 2) ? Pi/2.0 : (2.0*Pi)/fNtheta;
	double angle  = theta - dtheta;
	
	for (int i = 0; i < fNtheta; i++)
	{
		/* orientation */
		angle += dtheta;
	
		/* components */
		fAngles[i] = angle;
	}

	return fAngles;
}
