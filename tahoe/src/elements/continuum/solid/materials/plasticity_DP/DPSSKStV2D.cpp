/* $Id: DPSSKStV2D.cpp,v 1.9.4.1 2004-04-08 07:33:08 paklein Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
DPSSKStV2D::DPSSKStV2D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_StVenant_DP_2D"),
	DPSSKStV(in, support),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{

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
	return fModulus2D;
}

/* stress */
const dSymMatrixT& DPSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(DPSSKStV::s_ij());
	return fStress2D;
}

/* describe the parameters needed by the interface */
void DPSSKStV2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	DPSSKStV::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}
