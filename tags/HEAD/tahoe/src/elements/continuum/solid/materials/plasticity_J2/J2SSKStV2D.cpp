/* $Id: J2SSKStV2D.cpp,v 1.1.1.1 2001-01-29 08:20:30 paklein Exp $ */
/* created: paklein (06/18/1997)                                          */

#include "J2SSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

/* constructor */
J2SSKStV2D::J2SSKStV2D(ifstreamT& in, const ElasticT& element):
J2SSKStV(in, element),
	Material2DT(in, Material2DT::kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{

}

/* returns elastic strain (3D) */
const dSymMatrixT& J2SSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return J2SSKStV::ElasticStrain(fTotalStrain3D, element, ip);
}

/* print parameters */
void J2SSKStV2D::Print(ostream& out) const
{
	/* inherited */
	J2SSKStV::Print(out);
	Material2DT::Print(out);
}

/* moduli */
const dMatrixT& J2SSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(J2SSKStV::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

/* stress */
const dSymMatrixT& J2SSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(J2SSKStV::s_ij());
	fStress2D *= fThickness;
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double J2SSKStV2D::StrainEnergyDensity(void)
{
	return fThickness*J2SSKStV::StrainEnergyDensity();
}

/***********************************************************************
* Protected
***********************************************************************/

/* print name */
void J2SSKStV2D::PrintName(ostream& out) const
{
	/* inherited */
	J2SSKStV::PrintName(out);
	out << "    2D\n";
}
