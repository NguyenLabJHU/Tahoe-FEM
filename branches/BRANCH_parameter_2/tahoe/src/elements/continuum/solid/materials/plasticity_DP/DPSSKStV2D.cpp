/* $Id: DPSSKStV2D.cpp,v 1.8.32.1 2004-01-21 19:10:20 paklein Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
DPSSKStV2D::DPSSKStV2D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_StVenant_DP_2D"),
	DPSSKStV(in, support),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* initialization */
void DPSSKStV2D::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* returns elastic strain (3D) */
const dSymMatrixT& DPSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return DPSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

}

/* print parameters */
void DPSSKStV2D::Print(ostream& out) const
{
	/* inherited */
	DPSSKStV::Print(out);
	Material2DT::Print(out);
}

/* print name */
void DPSSKStV2D::PrintName(ostream& out) const
{
	/* inherited */
	DPSSKStV::PrintName(out);
	out << "    2D\n";
}

/* moduli */
const dMatrixT& DPSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(DPSSKStV::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& DPSSKStV2D::cdisc_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(DPSSKStV::cdisc_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}


/* stress */
const dSymMatrixT& DPSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(DPSSKStV::s_ij());
	fStress2D *= fThickness;  
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double DPSSKStV2D::StrainEnergyDensity(void)
{
	return fThickness*DPSSKStV::StrainEnergyDensity();
}
