/* $Id: SimoIso2D.cpp,v 1.4 2001-07-03 01:35:14 paklein Exp $ */
/* created: paklein (03/04/1997)                                          */
/* (2D <-> 3D) translator for the SimoIso3D.                              */

#include "SimoIso2D.h"
#include <math.h>
#include <iostream.h>

/* constructor */
SimoIso2D::SimoIso2D(ifstreamT& in, const FiniteStrainT& element):
	SimoIso3D(in, element),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fb_2D(2)
{
	fDensity *= fThickness;
}

/* moduli */
const dMatrixT& SimoIso2D::c_ijkl(void)
{
	/* b */
	Compute_b(fb_2D);
	
	/* Compute plane strain stretch */
	fb.ExpandFrom2D(fb_2D);
	fb(2,2) = 1.0; /* plane strain */

	/* compute b_bar */
	double J = fb_2D.Det();
	if (J <= 0.0) throw eBadJacobianDet;
	J = sqrt(J);
	fb_bar.SetToScaled(pow(J,-2.0/3.0), fb);

	/* 3D calculation */
	ComputeModuli(J, fb_bar, fModulus);

	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(fModulus);
	fModulus2D *= fThickness;

	return fModulus2D;
}
	
/* stresses */
const dSymMatrixT& SimoIso2D::s_ij(void)
{
	/* b */
	Compute_b(fb_2D);

	/* Compute plane strain stretch */
	fb.ExpandFrom2D(fb_2D);
	fb(2,2) = 1.0; //out-of-plane stretch
	
	/* compute b_bar */
	double J = fb_2D.Det();
	if (J <= 0.0) throw eBadJacobianDet;
	J = sqrt(J);
	fb_bar.SetToScaled(pow(J,-2.0/3.0), fb);

	/* 3D calculation */
	ComputeCauchy(J, fb_bar, fStress);

	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(fStress);
	fStress2D *= fThickness;

	return fStress2D;
}

/* strain energy density */
double SimoIso2D::StrainEnergyDensity(void)
{
	/* b */
	Compute_b(fb_2D);

	/* Compute plane strain stretch */
	fb.ExpandFrom2D(fb_2D);
	fb(2,2) = 1.0; //out-of-plane stretch
	
	/* compute b_bar */
	double J = fb_2D.Det();
	if (J <= 0.0) throw eBadJacobianDet;
	J = sqrt(J);
	fb_bar.SetToScaled(pow(J,-2.0/3.0), fb);

	return fThickness*ComputeEnergy(J, fb);
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
