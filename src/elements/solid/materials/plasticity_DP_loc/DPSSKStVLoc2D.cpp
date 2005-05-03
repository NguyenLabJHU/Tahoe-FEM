/* $Id: DPSSKStVLoc2D.cpp,v 1.8 2005-05-03 15:52:28 raregue Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStVLoc2D.h"
#include "SSEnhLocMatSupportT.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "DPSSLinHardLocT.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

using namespace Tahoe;

/* constructor */
DPSSKStVLoc2D::DPSSKStVLoc2D(void):
	ParameterInterfaceT("small_strain_StVenant_DP_Loc_2D"),
	fSSEnhLocMatSupport(NULL)
{

}

/* returns elastic strain (3D) */
const dSymMatrixT& DPSSKStVLoc2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return DPSSKStVLoc::ElasticStrain(fTotalStrain3D, element, ip);
}

/* tangent modulus */
const dMatrixT& DPSSKStVLoc2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(DPSSKStVLoc::c_ijkl());
//	fModulus2D *= fThickness;
	return fModulus2D;
}

/* elastic modulus */
const dMatrixT& DPSSKStVLoc2D::ce_ijkl(void)
{
	/* 3D -> 2D */
	fModulusElas2D.Rank4ReduceFrom3D(DPSSKStVLoc::ce_ijkl());
//	fModulus2D *= fThickness;
	return fModulusElas2D;
}

/* perfectly-plastic continuum modulus */
const dMatrixT& DPSSKStVLoc2D::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulusPerfPlas2D.Rank4ReduceFrom3D(DPSSKStVLoc::c_perfplas_ijkl());
//	fModulusPerfPlas2D *= fThickness;
	return fModulusPerfPlas2D;
}


/* stress */
const dSymMatrixT& DPSSKStVLoc2D::s_ij(void)
{
#ifdef ENHANCED_STRAIN_LOC_DEV	
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	int element_locflag = 0;
	if (element.IsAllocated()) 
	{
		element_locflag = fSSEnhLocMatSupport->ElementLocflag();
	}
	if ( element_locflag == 2 )
	{
		fStress2D = fSSEnhLocMatSupport->ElementStress(ip);
	}
	else
	{
		/* 3D -> 2D */
		fStress2D.ReduceFrom3D(DPSSKStVLoc::s_ij());
	}	
#else
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(DPSSKStVLoc::s_ij());
#endif
	return fStress2D;
}

/* describe the parameters needed by the interface */
void DPSSKStVLoc2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	DPSSKStVLoc::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void DPSSKStVLoc2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	DPSSKStVLoc::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusElas2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
	
	/* cast to small strain embedded discontinuity material pointer */
	fSSEnhLocMatSupport = TB_DYNAMIC_CAST(const SSEnhLocMatSupportT*, fSSMatSupport);
}
