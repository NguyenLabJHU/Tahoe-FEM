/* $Id: GradJ2SSNonlinHard2D.cpp,v 1.3 2002-11-14 17:06:29 paklein Exp $ */
#include "GradJ2SSNonlinHard2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
GradJ2SSNonlinHard2D::GradJ2SSNonlinHard2D(ifstreamT& in, const SSMatSupportT& support) :
  GradJ2SSNonlinHard(in, support),  
  Material2DT(in, Material2DT::kPlaneStrain),
  fStress2D(2),
  fModulus2D(dSymMatrixT::NumValues(2)),
  fTotalStrain3D(3)
{
	/* acccount for thickness */
	fDensity *= fThickness;
}

/* initialization */
void GradJ2SSNonlinHard2D::Initialize(void)
{
	/* inherited */
	GradJ2SSNonlinHard::Initialize();
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

/* print parameters */
void GradJ2SSNonlinHard2D::Print(ostream& out) const
{
	/* inherited */
	GradJ2SSNonlinHard::Print(out);
	Material2DT::Print(out);
}

/* moduli */
const dMatrixT& GradJ2SSNonlinHard2D::c_ijkl()
{
	/* 3D -> 2D */
        fModulus2D.Rank4ReduceFrom3D(GradJ2SSNonlinHard::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

/* stress */
const dSymMatrixT& GradJ2SSNonlinHard2D::s_ij()
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(GradJ2SSNonlinHard::s_ij());
	fStress2D *= fThickness;
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double GradJ2SSNonlinHard2D::StrainEnergyDensity(void)
{
	return fThickness*GradJ2SSNonlinHard::StrainEnergyDensity();
}

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void GradJ2SSNonlinHard2D::PrintName(ostream& out) const
{
  // inherited
  GradJ2SSNonlinHard::PrintName(out);

  // output model name
  out << "    Plane Strain\n";
}
