/* $Id: FossumSSIso2DT.cpp,v 1.5 2004-05-08 23:46:35 raregue Exp $ */
#include "FossumSSIso2DT.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
FossumSSIso2DT::FossumSSIso2DT(ifstreamT& in, const SSMatSupportT& support):
	FossumSSIsoT(in, support),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* initialization */
void FossumSSIso2DT::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* returns elastic strain (3D) */
const dSymMatrixT& FossumSSIso2DT::ElasticStrain(const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return FossumSSIsoT::ElasticStrain(fTotalStrain3D, element, ip);
}

/* print parameters */
void FossumSSIso2DT::Print(ostream& out) const
{
	/* inherited */
	FossumSSIsoT::Print(out);
	Material2DT::Print(out);
}

/* print name */
void FossumSSIso2DT::PrintName(ostream& out) const
{
	/* inherited */
	FossumSSIsoT::PrintName(out);
	out << " 2D\n";
}

/* moduli */
const dMatrixT& FossumSSIso2DT::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(FossumSSIsoT::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& FossumSSIso2DT::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(FossumSSIsoT::c_perfplas_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}


/* stress */
const dSymMatrixT& FossumSSIso2DT::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(FossumSSIsoT::s_ij());
	fStress2D *= fThickness;  
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double FossumSSIso2DT::StrainEnergyDensity(void)
{
	return fThickness*FossumSSIsoT::StrainEnergyDensity();
}
