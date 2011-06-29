/* $Id: SimoIso2D.cpp,v 1.11 2004-09-10 22:39:12 paklein Exp $ */
/* created: paklein (03/04/1997) */
#include "SimoIso2D.h"
#include <math.h>

using namespace Tahoe;

/* constructor */
SimoIso2D::SimoIso2D(void):
	ParameterInterfaceT("Simo_isotropic_2D")
{
	/* set default value */
	fConstraint = kPlaneStrain;
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
			ExceptionT::GeneralFail("SimoIso2D::InitStep", "expecting isotropic (F_thermal)^-1:");
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

	return ComputeEnergy(J, fb);
}

/* accept parameter list */
void SimoIso2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SimoIso3D::TakeParameterList(list);
	
	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fb_2D.Dimension(2);	
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
