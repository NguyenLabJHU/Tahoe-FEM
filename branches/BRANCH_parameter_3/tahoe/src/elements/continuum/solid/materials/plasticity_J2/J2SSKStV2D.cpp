/* $Id: J2SSKStV2D.cpp,v 1.4.48.2 2004-06-08 22:27:33 paklein Exp $ */
/* created: paklein (06/18/1997) */
#include "J2SSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
J2SSKStV2D::J2SSKStV2D(ifstreamT& in, const SSMatSupportT& support):
	ParameterInterfaceT("small_strain_J2_StVenant_2D"),
	J2SSKStV(in, support),
	fStress2D(2),
	fModulus2D(dSymMatrixT::NumValues(2)),
	fTotalStrain3D(3)
{

}

J2SSKStV2D::J2SSKStV2D(void):
	ParameterInterfaceT("small_strain_J2_StVenant_2D")
{

}

/* initialization */
void J2SSKStV2D::Initialize(void)
{
	/* inherited */
	HookeanMatT::Initialize();
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

/* moduli */
const dMatrixT& J2SSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(J2SSKStV::c_ijkl());
	return fModulus2D;
}

/* stress */
const dSymMatrixT& J2SSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(J2SSKStV::s_ij());
	return fStress2D;
}

/* describe the parameters needed by the interface */
void J2SSKStV2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	J2SSKStV::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void J2SSKStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	J2SSKStV::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
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
