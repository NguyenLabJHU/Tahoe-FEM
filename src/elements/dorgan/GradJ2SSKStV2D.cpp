/* $Id: GradJ2SSKStV2D.cpp,v 1.1 2004-09-02 18:25:04 rdorgan Exp $ */
#include "GradJ2SSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"

using namespace Tahoe;

/* constructor */
GradJ2SSKStV2D::GradJ2SSKStV2D(void):
	ParameterInterfaceT("grad_small_strain_StVenant_J2_2D")
{

}

/* returns elastic strain (3D) */
const dSymMatrixT& GradJ2SSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain,
	const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return GradJ2SSKStV::ElasticStrain(fTotalStrain3D, element, ip);
}

/* moduli */
const dMatrixT& GradJ2SSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(GradJ2SSKStV::c_ijkl());
	return fModulus2D;
}

/* off diagonal moduli for Kar */
const dMatrixT& GradJ2SSKStV2D::odm_bh_ij(void)
{
	fOffDiagonalModulus_bh_2D = 0;

	/* 3D -> 2D */
	ReduceOffDiagonalModulus_bh(GradJ2SSKStV::odm_bh_ij(), fOffDiagonalModulus_bh_2D);

	/* 3D -> 2D */
	//	fModulus2D.Rank4ReduceFrom3D(GradJ2SSKStV::c_ijkl());
	//	fStress2D.ReduceFrom3D(GradJ2SSKStV::n_ij());

	//	fModulus2D.Multx(fStress2D,fTensorTemp1);
	//	fOffDiagonalModulus_bh_2D.AddScaled(-1.0, fTensorTemp1);
	
	return fOffDiagonalModulus_bh_2D;
}

/* off diagonal moduli for Kar */
const dMatrixT& GradJ2SSKStV2D::odm_hb_ij(void)
{
	fOffDiagonalModulus_hb_2D.Transpose(GradJ2SSKStV2D::odm_bh_ij());
	return fOffDiagonalModulus_hb_2D;
}

/* moduli for gradient term in K_hp */
const dMatrixT& GradJ2SSKStV2D::gm_hp(void)
{
	/* 3D -> 2D */
	ReduceGradientModuli_hp(GradJ2SSKStV::gm_hp(), fGradientModulus_hp_2D);
	return fGradientModulus_hp_2D;
}

/* stress */
const dSymMatrixT& GradJ2SSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(GradJ2SSKStV::s_ij());
	return fStress2D;
}

/* describe the parameters needed by the interface */
void GradJ2SSKStV2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	GradJ2SSKStV::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void GradJ2SSKStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	GradJ2SSKStV::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);

	fOffDiagonalModulus_bh_2D.Dimension(dSymMatrixT::NumValues(2),1);
	fOffDiagonalModulus_hb_2D.Dimension(1,dSymMatrixT::NumValues(2));
	fGradientModulus_hp_2D.Dimension(1,2);
	fTensorTemp1.Dimension(dSymMatrixT::NumValues(2),1);
}
