/* $Id: LatLongPtsT.cpp,v 1.1.1.1.10.1 2002-06-27 18:03:24 cjkimme Exp $ */
/* created: paklein (10/31/1997)                                          */
/* Base class for spherical point generators.                             */

#include "LatLongPtsT.h"

#include <math.h>
#include <iostream.h>

#include "Constants.h"
#include "ExceptionCodes.h"
#include "fstreamT.h"


using namespace Tahoe;

const double Pi = acos(-1.0);

/*
* Constructor
*/
LatLongPtsT::LatLongPtsT(ifstreamT& in)
{
	/* number of integration points */
	in >> fNphi >> fNtheta;
	if (fNphi < 1 || fNtheta < 1) throw eBadInputValue;
	
	int Ntot = fNtheta*fNphi + 2;  //2 extra for the z-axis caps
	fPoints.Allocate(Ntot,3);
	fJacobians.Allocate(Ntot);
}

/*
* Print parameters.
*/
void LatLongPtsT::Print(ostream& out) const
{
	/* number of integration points */
	out << " Number of phi sampling points . . . . . . . . . = " << fNphi   << '\n';
	out << " Number of theta sampling points . . . . . . . . = " << fNtheta << '\n';
}

void LatLongPtsT::PrintName(ostream& out) const
{
	out << "    Lattitude and Longitude lines\n";
}

/*
* Generate sphere points:
*
*   theta = angle about z from x
*   phi   = angle about x from z
*
* The final orientation is generated by applied the
* theta and phi rotations in succession about the local
* axes.
*
*/
const dArray2DT& LatLongPtsT::SpherePoints(double phi_tr, double theta_tr)
{
	/* fill in angle tables */
	double dphi    = (2.0*Pi)/fNphi;	
	double dtheta  = Pi/fNtheta;
	double jfactor = dphi*dtheta;
	double theta  =-dtheta/2.0;

	int    thetacount = fNtheta - 1;
	double jsum = 0.0;
	
	dArrayT	xsi;
	int Ntot = fJacobians.Length();
	for (int i = 0; i < Ntot-2; i++) //skip caps
	{
		double phi;
		if (++thetacount == fNtheta)
		{
			phi    = 0.0;
			theta += dtheta;
			thetacount = 0;			
		}
		else
			phi += dphi;

		/* jacobian */
		fJacobians[i] = jfactor*sin(theta);
		jsum += fJacobians[i];

		/* direction cosines */
		fPoints.RowAlias(i,xsi);
		xsi[0] = sin(theta)*cos(phi);
		xsi[1] = sin(theta)*sin(phi);
		xsi[2] = cos(theta);
	}

	double xsi3 = 1.0;
	for (int j = Ntot-2; j < Ntot; j++) //caps
	{
		/* jacobian */
		//fjacobian[j] = 2.0*Pi*sin(dtheta/2.0);  //actual cap area
		fJacobians[j] = (4.0*Pi - jsum)/2.0;      //correct total area

		/* direction cosines */
		fPoints.RowAlias(j,xsi);
		xsi[0] = 0.0;
		xsi[1] = 0.0;
		xsi[2] = xsi3;
	}

	/* reorient points */
	TransformPoints(phi_tr,theta_tr);
	
	return fPoints;
}
