/* $Id: DPSSKStVLoc2D.cpp,v 1.5 2005-02-16 17:26:24 raregue Exp $ */
/* created: myip (06/01/1999) */
#include "DPSSKStVLoc2D.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "DPSSLinHardLocT.h"

using namespace Tahoe;

/* constructor */
DPSSKStVLoc2D::DPSSKStVLoc2D(void):
	ParameterInterfaceT("small_strain_StVenant_DP_Loc_2D")
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
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(DPSSKStVLoc::s_ij());
//	fStress2D *= fThickness;  
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
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
}
