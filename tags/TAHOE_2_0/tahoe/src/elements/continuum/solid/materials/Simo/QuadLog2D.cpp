/* $Id: QuadLog2D.cpp,v 1.6 2004-07-15 08:27:35 paklein Exp $ */
/* created: paklein (06/28/1997) */
#include "QuadLog2D.h"
#include <math.h>
#include <iostream.h>

using namespace Tahoe;

/* constructor */
QuadLog2D::QuadLog2D(void): ParameterInterfaceT("quad_log_2D") { }

/* modulus */
const dMatrixT& QuadLog2D::c_ijkl(void)
{
	/* deformation */
	Compute_b(fb_2D);

	/* Compute plane strain stretch */
	fb.ExpandFrom2D(fb_2D);
	fb(2,2) = 1.0; /* plane strain */
	
	/* 3D calculation */
	ComputeModuli(fb, fModulus);

	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(fModulus);

	return fModulus2D;
}
	
/* stresses */
const dSymMatrixT& QuadLog2D::s_ij(void)
{
	/* deformation */
	Compute_b(fb_2D);

	/* Compute plane strain stretch */
	fb.ExpandFrom2D(fb_2D);
	fb(2,2) = 1.0; /* plane strain */
	
	/* 3D calculation */
	ComputeCauchy(fb, fStress);

	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(fStress);

	return fStress2D;
}

/* strain energy density */
double QuadLog2D::StrainEnergyDensity(void)
{
	/* deformation */
	Compute_b(fb_2D);

	/* principal values */
	fb_2D.PrincipalValues(fEigs);
	fEigs[2] = 1.0; /* plane strain */
	
	/* logarithmic stretches */
	LogStretches(fEigs);

	return ComputeEnergy(floge);
}

/* describe the parameters needed by the interface */
void QuadLog2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	QuadLog3D::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void QuadLog2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	QuadLog3D::TakeParameterList(list);

	fb_2D.Dimension(2);
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
}
