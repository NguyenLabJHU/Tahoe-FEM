/* $Id: J2SSKStV2D.cpp,v 1.6 2004-09-10 22:39:32 paklein Exp $ */
/* created: paklein (06/18/1997) */
#include "J2SSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
J2SSKStV2D::J2SSKStV2D(void):
	ParameterInterfaceT("small_strain_StVenant_J2_2D")
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

/* returns elastic strain (3D) */
const dSymMatrixT& J2SSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int nip, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return J2SSKStV::ElasticStrain(fTotalStrain3D, element, nip, ip);
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
