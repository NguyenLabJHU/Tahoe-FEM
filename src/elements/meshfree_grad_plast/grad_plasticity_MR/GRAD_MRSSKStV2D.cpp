/* created: Majid T. Manzari (04/16/2003) */
#include "GRAD_MRSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
GRAD_MRSSKStV2D::GRAD_MRSSKStV2D(ifstreamT& in, const SSMatSupportT& support):
	GRAD_MRSSKStV(in, support),
	Material2DT(in, kPlaneStrain),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fYieldFunction2D(0.0),
	fTotalStrain3D(3)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* initialization */
void GRAD_MRSSKStV2D::Initialize(void)
{
ExceptionT::GeneralFail("GRAD_MRSSKStV2D::Initialize", "out of date");
#if 0
	/* inherited */
	HookeanMatT::Initialize();
#endif
}

/* returns 3D total strain (3D) */
const dSymMatrixT& GRAD_MRSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip) //del2_totalstrain??
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	/*return fTotalStrain3D;*/
	return GRAD_MRSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

}

/* returns 3D  gradient of total strain (3D) */
const dSymMatrixT& GRAD_MRSSKStV2D::GradElasticStrain(const dSymMatrixT& del2_totalstrain, 
	const ElementCardT& element, int ip) //del2_totalstrain??
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(del2_totalstrain);

	/* inherited */
	/*return fTotalStrain3D;*/
	return GRAD_MRSSKStV::GradElasticStrain(fTotalStrain3D, element, ip);

}

/* print parameters */
void GRAD_MRSSKStV2D::Print(ostream& out) const
{
	/* inherited */
	GRAD_MRSSKStV::Print(out);
	Material2DT::Print(out);
}

/* print name */
void GRAD_MRSSKStV2D::PrintName(ostream& out) const
{
	/* inherited */
	GRAD_MRSSKStV::PrintName(out);
	out << "    2D\n";
}

/* moduli */
const dMatrixT& GRAD_MRSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::c_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& GRAD_MRSSKStV2D::cdisc_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::cdisc_ijkl());
	fModulus2D *= fThickness;
	return fModulus2D;
}


/* stress */
const dSymMatrixT& GRAD_MRSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(GRAD_MRSSKStV::s_ij());
	fStress2D *= fThickness;  
	return fStress2D;
}

/* yield function */
const double& GRAD_MRSSKStV2D::Yield_Function(void)
{
	
	fYieldFunction2D = GRAD_MRSSKStV::YieldF();
	return fYieldFunction2D;
}

/* returns the strain energy density for the specified strain */
double GRAD_MRSSKStV2D::StrainEnergyDensity(void)
{
	return fThickness*GRAD_MRSSKStV::StrainEnergyDensity();
}