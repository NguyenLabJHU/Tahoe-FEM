/* $Id: SimoIso2D.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (03/04/1997)                                          */
/* (2D <-> 3D) translator for the SimoIso3D.                              */

#include "SimoIso2D.h"
#include <math.h>
#include <iostream.h>

/* constructor */
SimoIso2D::SimoIso2D(ifstreamT& in, const ElasticT& element):
	SimoIso3D(in, element),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fb_3D(3)
{
	fDensity *= fThickness;
}

/* moduli */
const dMatrixT& SimoIso2D::c_ijkl(void)
{
	/* Compute plane strain stretch */
	fb_3D.ExpandFrom2D(b());
	fb_3D(2,2) = 1.0; //out-of-plane stretch
	
	/* 3D calculation */
	double J = sqrt(fb_3D.Det());
	fb_3D *= pow(J,-2.0/3.0);

	ComputeModuli(J, fb_3D, fModulus);

	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(fModulus);
	fModulus2D *= fThickness;

	return fModulus2D;
}
	
/* stresses */
const dSymMatrixT& SimoIso2D::s_ij(void)
{
	/* Compute plane strain stretch */
	fb_3D.ExpandFrom2D(b());
	fb_3D(2,2) = 1.0; //out-of-plane stretch
	
	/* 3D calculation */
	double J = sqrt(fb_3D.Det());
	fb_3D *= pow(J,-2.0/3.0);

	ComputeCauchy(J, fb_3D, fStress);

	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(fStress);
	fStress2D *= fThickness;

	return fStress2D;
}

/* strain energy density */
double SimoIso2D::StrainEnergyDensity(void)
{
	/* Compute plane strain stretch */
	fb_3D.ExpandFrom2D(b());
	fb_3D(2,2) = 1.0; //out-of-plane stretch
	
	/* 3D calculation */
	double J = sqrt(fb_3D.Det());
	fb_3D *= pow(J,-2.0/3.0);

	return fThickness*ComputeEnergy(J, fb_3D);
}

/* print parameters */
void SimoIso2D::Print(ostream& out) const
{
	/* inherited */
	SimoIso3D::Print(out);
	Material2DT::Print(out);
}

/*************************************************************************
* Protected
*************************************************************************/

/* print name */
void SimoIso2D::PrintName(ostream& out) const
{
	/* inherited */
	SimoIso3D::PrintName(out);

	out << "    Plane Strain\n";
}
