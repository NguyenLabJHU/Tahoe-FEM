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
	fTotalStrain3D(3)
{
	/* account for thickness */
	fDensity *= fThickness;
}

/* initialization */
void GRAD_MRSSKStV2D::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
}

/* returns 3D total strain (3D) */
const dSymMatrixT& GRAD_MRSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip) //gradtotalstrain??
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	/*return fTotalStrain3D;*/
	return GRAD_MRSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

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

//include yield function here
const double GRAD_MRSSKStV2D::yield_function(void)
{
	
	fYieldFunction2D = fYieldFunction;
	return fYieldFunction2D;
}

/* returns the strain energy density for the specified strain */
double GRAD_MRSSKStV2D::StrainEnergyDensity(void)
{
	return fThickness*GRAD_MRSSKStV::StrainEnergyDensity();
}