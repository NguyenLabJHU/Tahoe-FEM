/* created: Majid T. Manzari (04/16/2003) */
#include "MRSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
MRSSKStV2D::MRSSKStV2D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("MRSSKStV2D"),
	MRSSKStV(in, support),
//	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{
	/* account for thickness */
//	fDensity *= fThickness;
}

/* initialization */
void MRSSKStV2D::Initialize(void)
{
ExceptionT::GeneralFail("MRSSKStV2D::Initialize", "out of date");
#if 0
	/* inherited */
	HookeanMatT::Initialize();
#endif
}

/* returns 3D total strain (3D) */
const dSymMatrixT& MRSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	/*return fTotalStrain3D;*/
	return MRSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

}

/* moduli */
const dMatrixT& MRSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(MRSSKStV::c_ijkl());
//	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& MRSSKStV2D::cdisc_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(MRSSKStV::cdisc_ijkl());
//	fModulus2D *= fThickness;
	return fModulus2D;
}


/* stress */
const dSymMatrixT& MRSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(MRSSKStV::s_ij());
//	fStress2D *= fThickness;  
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double MRSSKStV2D::StrainEnergyDensity(void)
{
	return MRSSKStV::StrainEnergyDensity();
}