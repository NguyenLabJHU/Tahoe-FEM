/* $Id: MRSSKStV2D.cpp,v 1.6 2010-07-21 19:58:20 regueiro Exp $ */
/* created: Majid T. Manzari (04/16/2003) */
#include "MRSSKStV2D.h"
#include "SSEnhLocMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "MRSSNLHardT.h"

#include "DevelopmentElementsConfig.h"

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

/* elastic modulus */
const dMatrixT& MRSSKStV2D::ce_ijkl(void)
{
	/* 3D -> 2D */
	fModulusElas2D.Rank4ReduceFrom3D(MRSSKStV::ce_ijkl());
	return fModulusElas2D;
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
#ifdef ENHANCED_STRAIN_LOC_DEV	
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	int elem = CurrElementNumber();
	element_locflag = 0;
	if (element.IsAllocated()) element_locflag = fSSEnhLocMatSupport->ElementLocflag(elem);
	if ( element_locflag == 2 )
	{
		fStress2D = fSSEnhLocMatSupport->ElementStress(elem,ip);
	}
	else
	{
		/* 3D -> 2D */
		fStress2D.ReduceFrom3D(MRSSKStV::s_ij());
	}	
#else
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(MRSSKStV::s_ij());
#endif
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
	fModulusElas2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
}