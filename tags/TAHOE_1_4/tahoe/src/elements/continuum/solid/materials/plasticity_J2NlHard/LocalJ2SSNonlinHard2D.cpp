/* $Id: LocalJ2SSNonlinHard2D.cpp,v 1.3 2002-11-14 17:06:29 paklein Exp $ */
#include "LocalJ2SSNonlinHard2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
LocalJ2SSNonlinHard2D::LocalJ2SSNonlinHard2D(ifstreamT& in, const SSMatSupportT& support) :
  LocalJ2SSNonlinHard(in, support),  
  Material2DT(in, Material2DT::kPlaneStrain),
  fStress2D(2),
  fModulus2D(dSymMatrixT::NumValues(2)),
  fTotalStrain3D(3)
{
	/* acccount for thickness */
	fDensity *= fThickness;
}

/* initialization */
void LocalJ2SSNonlinHard2D::Initialize(void)
{
	/* inherited */
	LocalJ2SSNonlinHard::Initialize();
}

/* returns elastic strain (3D) */
const dSymMatrixT& LocalJ2SSNonlinHard2D::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return LocalJ2SSNonlinHard::ElasticStrain(fTotalStrain3D, element, ip);
}

/* print parameters */
void LocalJ2SSNonlinHard2D::Print(ostream& out) const
{
	/* inherited */
	LocalJ2SSNonlinHard::Print(out);
	Material2DT::Print(out);
}

/* moduli */
const dMatrixT& LocalJ2SSNonlinHard2D::c_ijkl()
{
	/* 3D -> 2D */
        fModulus2D.Rank4ReduceFrom3D(LocalJ2SSNonlinHard::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

/* stress */
const dSymMatrixT& LocalJ2SSNonlinHard2D::s_ij()
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(LocalJ2SSNonlinHard::s_ij());
	fStress2D *= fThickness;
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double LocalJ2SSNonlinHard2D::StrainEnergyDensity(void)
{
	return fThickness*LocalJ2SSNonlinHard::StrainEnergyDensity();
}

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void LocalJ2SSNonlinHard2D::PrintName(ostream& out) const
{
  // inherited
  LocalJ2SSNonlinHard::PrintName(out);

  // output model name
  out << "    Plane Strain\n";
}
