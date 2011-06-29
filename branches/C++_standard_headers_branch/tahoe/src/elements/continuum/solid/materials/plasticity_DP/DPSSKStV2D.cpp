/* $Id: DPSSKStV2D.cpp,v 1.12 2005-02-25 18:41:18 cfoster Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "DPSSLinHardT.h"

using namespace Tahoe;

/* constructor */
DPSSKStV2D::DPSSKStV2D(void):
	ParameterInterfaceT("small_strain_StVenant_DP_2D")
{
	/* reset default value */
	fConstraint = kPlaneStrain;
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

/* moduli */
const dMatrixT& DPSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(DPSSKStV::c_ijkl());
	return fModulus2D;
}

/* moduli */
const dMatrixT& DPSSKStV2D::ce_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(DPSSKStV::ce_ijkl());
	return fModulus2D;
}

/* stress */
const dSymMatrixT& DPSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(DPSSKStV::s_ij());
	return fStress2D;
}

/* accept parameter list */
void DPSSKStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	DPSSKStV::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
}
