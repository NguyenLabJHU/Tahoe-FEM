/* $Id: LatLongPtsT.cpp,v 1.6 2004-07-15 08:28:09 paklein Exp $ */
/* created: paklein (10/31/1997) */
#include "LatLongPtsT.h"

#include <math.h>
#include "ExceptionT.h"

using namespace Tahoe;

const double Pi = acos(-1.0);

/* constructor */
LatLongPtsT::LatLongPtsT(int n_phi, int n_theta):
	fNphi(n_phi),
	fNtheta(n_theta)
{
	if (fNphi < 1 || fNtheta < 1) ExceptionT::BadInputValue("LatLongPtsT::LatLongPtsT");
	
	int Ntot = fNtheta*fNphi + 2;  //2 extra for the z-axis caps
	fPoints.Dimension(Ntot,3);
	fJacobians.Dimension(Ntot);
}

/* generate sphere points */
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
