/* $Id: DPSSKStVLoc2D.cpp,v 1.2 2004-06-09 17:27:39 raregue Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStVLoc2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
DPSSKStVLoc2D::DPSSKStVLoc2D(ifstreamT& in, const SSMatSupportT& support):
	DPSSKStVLoc(in, support),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fModulusPerfPlas2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* initialization */
void DPSSKStVLoc2D::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* returns elastic strain (3D) */
const dSymMatrixT& DPSSKStVLoc2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return DPSSKStVLoc::ElasticStrain(fTotalStrain3D, element, ip);

}

/* print parameters */
void DPSSKStVLoc2D::Print(ostream& out) const
{
	/* inherited */
	DPSSKStVLoc::Print(out);
	Material2DT::Print(out);
}

/* print name */
void DPSSKStVLoc2D::PrintName(ostream& out) const
{
	/* inherited */
	DPSSKStVLoc::PrintName(out);
	out << "    2D\n";
}

/* moduli */
const dMatrixT& DPSSKStVLoc2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(DPSSKStVLoc::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& DPSSKStVLoc2D::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulusPerfPlas2D.Rank4ReduceFrom3D(DPSSKStVLoc::c_perfplas_ijkl());
	fModulusPerfPlas2D *= fThickness;
	return fModulusPerfPlas2D;
}


/* stress */
const dSymMatrixT& DPSSKStVLoc2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(DPSSKStVLoc::s_ij());
	fStress2D *= fThickness;  
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double DPSSKStVLoc2D::StrainEnergyDensity(void)
{
	return fThickness*DPSSKStVLoc::StrainEnergyDensity();
}
