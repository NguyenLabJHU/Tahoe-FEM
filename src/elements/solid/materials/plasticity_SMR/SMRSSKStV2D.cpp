/* $Id: SMRSSKStV2D.cpp,v 1.1 2006-07-27 13:20:08 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003) */
#include "SMRSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "SMRSSNLHardT.h"

using namespace Tahoe;

/* constructor */
SMRSSKStV2D::SMRSSKStV2D(void):
	ParameterInterfaceT("small_strain_StVenant_SMR_2D")
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

/* returns 3D total strain (3D) */
const dSymMatrixT& SMRSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return SMRSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

}

/* moduli */
const dMatrixT& SMRSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(SMRSSKStV::c_ijkl());
	return fModulus2D;
}

const dMatrixT& SMRSSKStV2D::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(SMRSSKStV::c_perfplas_ijkl());
	return fModulus2D;
}


/* stress */
const dSymMatrixT& SMRSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(SMRSSKStV::s_ij());  
	return fStress2D;
}

/* accept parameter list */
void SMRSSKStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SMRSSKStV::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
}