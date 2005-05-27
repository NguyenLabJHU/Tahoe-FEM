/* $Id: GradJ2SSNonlinHard2D.cpp,v 1.5 2004-09-10 22:39:37 paklein Exp $ */
#include "GradJ2SSNonlinHard2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
GradJ2SSNonlinHard2D::GradJ2SSNonlinHard2D(ifstreamT& in, const SSMatSupportT& support) :
	ParameterInterfaceT("small_strain_J2_nonlocal_2D"),
  GradJ2SSNonlinHard(in, support),  
  fStress2D(2),
  fModulus2D(dSymMatrixT::NumValues(2)),
  fTotalStrain3D(3)
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

/* returns elastic strain (3D) */
const dSymMatrixT& GradJ2SSNonlinHard2D::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return GradJ2SSNonlinHard::ElasticStrain(fTotalStrain3D, element, ip);
}

/* moduli */
const dMatrixT& GradJ2SSNonlinHard2D::c_ijkl()
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(GradJ2SSNonlinHard::c_ijkl());
	return fModulus2D;
}

/* stress */
const dSymMatrixT& GradJ2SSNonlinHard2D::s_ij()
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(GradJ2SSNonlinHard::s_ij());
	return fStress2D;
}
