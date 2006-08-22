/* $Id: MRSSKStV2D.cpp,v 1.5 2006-08-22 14:39:17 kyonten Exp $ */
/* created: Majid T. Manzari (04/16/2003) */
#include "MRSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "MRSSNLHardT.h"

using namespace Tahoe;

/* constructor */
MRSSKStV2D::MRSSKStV2D(void):
	ParameterInterfaceT("small_strain_StVenant_MR_2D")
{
	/* reset default value */
	fConstraint = kPlaneStrain;
}

/* returns 3D total strain (3D) */
const dSymMatrixT& MRSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return MRSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

}

/* moduli */
const dMatrixT& MRSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(MRSSKStV::c_ijkl());
	return fModulus2D;
}

const dMatrixT& MRSSKStV2D::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(MRSSKStV::c_perfplas_ijkl());
	return fModulus2D;
}


/* stress */
const dSymMatrixT& MRSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(MRSSKStV::s_ij());  
	return fStress2D;
}

/* accept parameter list */
void MRSSKStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	MRSSKStV::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
}