/* $Id: SimoIso2D.cpp,v 1.9 2003-01-29 07:34:48 paklein Exp $ */
/* created: paklein (03/04/1997) */
#include "SimoIso2D.h"
#include <math.h>
#include <iostream.h>

using namespace Tahoe;

/* constructor */
SimoIso2D::SimoIso2D(ifstreamT& in, const FSMatSupportT& support):
	SimoIso3D(in, support),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fb_2D(2)
{
	fDensity *= fThickness;
}

/* initialize step */
void SimoIso2D::InitStep(void)
{
	/* inherited */
	SimoIso3D::InitStep();

	/* check (inverse) thermal dilatation */
	const dMatrixT& F_therm_inv = F_thermal_inverse();
	if (HasThermalStrain())
	{
		/* inverse thermal dilatation */
		const dMatrixT& F_therm_inv = F_thermal_inverse();

		if (fabs(F_therm_inv(0,0) - F_therm_inv(1,1)) > kSmall ||
		    fabs(F_therm_inv(1,0)) > kSmall ||
		    fabs(F_therm_inv(0,1)) > kSmall)
		{
			cout << "\n SimoIso2D::InitStep: expecting isotropic (F_thermal)^-1:\n"
			     << F_therm_inv << endl;
			throw ExceptionT::kGeneralFail;
		}
	}
}

/* moduli */
const dMatrixT& SimoIso2D::c_ijkl(void)
{
	/* compute 3D stretch tensor */
	Compute_b_3D(fb);

	/* compute b_bar */
	double J = fb.Det();
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
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
	/* compute 3D stretch tensor */
	Compute_b_3D(fb);
	
	/* compute b_bar */
	double J = fb.Det();
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
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
	/* compute 3D stretch tensor */
	Compute_b_3D(fb);
	
	/* compute b_bar */
	double J = fb.Det();
	if (J <= 0.0) throw ExceptionT::kBadJacobianDet;
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

/*************************************************************************
* Private
*************************************************************************/

/** compute 3D stretch tensor \b b from the 2D deformation state. 
 * \todo Make this a FSSolidMatT function? */
void SimoIso2D::Compute_b_3D(dSymMatrixT& b_3D)
{
	/* get mechanical part of the deformation gradient */
	const dMatrixT& F_mech = F_mechanical();

	/* b */
	Compute_b(F_mech, fb_2D);
	
	/* Compute plane strain stretch */
	b_3D.ExpandFrom2D(fb_2D);
	if (HasThermalStrain()) /* assuming isotropic thermal strain */
	{
		double F_inv = (F_thermal_inverse())(0,0);
		b_3D(2,2) = F_inv*F_inv; 
	}
	else
		b_3D(2,2) = 1.0; /* plane strain */
}
