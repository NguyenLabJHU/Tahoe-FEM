/* created: Karma Yonten (03/04/2004)                   
   MR version modified to incorporate gradient plasticity 
   theory.
*/
#include "GRAD_MRSSKStV2D.h"
#include "ElementCardT.h"
#include "StringT.h"
#include "GRAD_MRSSNLHardT.h"

using namespace Tahoe;

/* constructor */
GRAD_MRSSKStV2D::GRAD_MRSSKStV2D(void):
	ParameterInterfaceT("small_strain_StVenant_GRAD_MR_2D")
{
	/* account for thickness */
//	fDensity *= fThickness;
}

/* returns 3D total strain (3D) */
const dSymMatrixT& GRAD_MRSSKStV2D::ElasticStrain(const dSymMatrixT& totalstrain, 
	const ElementCardT& element, int ip) 
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	/*return fTotalStrain3D;*/
	return GRAD_MRSSKStV::ElasticStrain(fTotalStrain3D, element, ip);

}

/* returns 3D  gradient of total strain (3D) */
const dSymMatrixT& GRAD_MRSSKStV2D::LapElasticStrain(const dSymMatrixT& laptotalstrain, 
	const ElementCardT& element, int ip) //lap_totalstrain??
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(laptotalstrain);

	/* inherited */
	/*return fTotalStrain3D;*/
	return GRAD_MRSSKStV::LapElasticStrain(fTotalStrain3D, element, ip);

}

/* moduli */
const dMatrixT& GRAD_MRSSKStV2D::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::c_ijkl());
//	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::c_perfplas_ijkl());
//	fModulus2D *= fThickness;
	return fModulus2D;
}


/* stress */
const dSymMatrixT& GRAD_MRSSKStV2D::s_ij(void)
{
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(GRAD_MRSSKStV::s_ij());
//	fStress2D *= fThickness;  
	return fStress2D;
}

/* yield function */
const double& GRAD_MRSSKStV2D::Yield_Function(void)
{
	
	fYieldFunction2D = GRAD_MRSSKStV::YieldF();
	return fYieldFunction2D;
}

/* describe the parameters needed by the interface */
void GRAD_MRSSKStV2D::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	GRAD_MRSSKStV::DefineParameters(list);
	
	/* 2D option must be plain stress */
	ParameterT& constraint = list.GetParameter("constraint_2D");
	constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void GRAD_MRSSKStV2D::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	GRAD_MRSSKStV::TakeParameterList(list);

	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2));
	fTotalStrain3D.Dimension(3);
	// fYieldFunction2D(0.0);  // scalar
	
}