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
	ParameterInterfaceT("small_strain_StVenant_MR_grad_2D")
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
	//fModulus2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::c_ijkl());
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

const dMatrixT& GRAD_MRSSKStV2D::c_UU1_ijkl(void)
{
	/* 3D -> 2D */
	fModulusUU1_2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::c_UU1_ijkl());
	return fModulusUU1_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_UU2_ijkl(void)
{
	/* 3D -> 2D */
	fModulusUU2_2D.Rank4ReduceFrom3D(GRAD_MRSSKStV::c_UU2_ijkl());
	return fModulusUU2_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_ULam1_ij(void)
{
	 fModulusULam1_2D = 0.0;
	/* 3D -> 2D */
	ReduceOffDiagonalModulus(GRAD_MRSSKStV::c_ULam1_ij(), fModulusULam1_2D);
	return fModulusULam1_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_ULam2_ij(void)
{
	fModulusULam2_2D = 0.0;
	/* 3D -> 2D */
	ReduceOffDiagonalModulus(GRAD_MRSSKStV::c_ULam2_ij(), fModulusULam2_2D);
	return fModulusULam2_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_LamU1_ij(void)
{
	fModulusLamU1_2D = 0.0;
	/* 3D -> 2D */
	fTemp1.Transpose(GRAD_MRSSKStV::c_LamU1_ij());
	ReduceOffDiagonalModulus(fTemp1, fTemp2);
	fModulusLamU1_2D.Transpose(fTemp2);
	return fModulusLamU1_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_LamU2_ij(void)
{
	fModulusLamU2_2D = 0.0;
	/* 3D -> 2D */
	fTemp3.Transpose(GRAD_MRSSKStV::c_LamU2_ij());
	ReduceOffDiagonalModulus(fTemp3, fTemp4);
	fModulusLamU2_2D.Transpose(fTemp4);
	return fModulusLamU2_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_LamLam1(void)
{
	/* 3D -> 2D */
	fModulusLamLam1_2D = GRAD_MRSSKStV::c_LamLam1();
	return fModulusLamLam1_2D;
}

const dMatrixT& GRAD_MRSSKStV2D::c_LamLam2(void)
{
	/* 3D -> 2D */
	fModulusLamLam2_2D = GRAD_MRSSKStV::c_LamLam2();
	return fModulusLamLam2_2D;
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
	fModulusUU1_2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusUU2_2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusULam1_2D.Dimension(dSymMatrixT::NumValues(2),1);
	fModulusULam2_2D.Dimension(dSymMatrixT::NumValues(2),1);
	fModulusLamU1_2D.Dimension(1,dSymMatrixT::NumValues(2));
	fModulusLamU2_2D.Dimension(1,dSymMatrixT::NumValues(2));
	fModulusLamLam1_2D.Dimension(1,1);
	fModulusLamLam2_2D.Dimension(1,1);
	fTotalStrain3D.Dimension(3);
	fTemp1.Dimension(dSymMatrixT::NumValues(2),1);
	fTemp2.Dimension(dSymMatrixT::NumValues(2),1);
	fTemp3.Dimension(dSymMatrixT::NumValues(2),1);
	fTemp4.Dimension(dSymMatrixT::NumValues(2),1);
	// fYieldFunction2D(0.0);  // scalar
	
}