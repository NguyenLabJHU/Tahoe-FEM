#include "GeoModelSS2DT.h"

#include "SSMatSupportT.h"
//#include "ElementCardT.h"
//#include "StringT.h"


using namespace Tahoe;

/* constructor */
GeoModelSS2DT::GeoModelSS2DT(void):
	ParameterInterfaceT("small_strain_GeoModel_2D"),
	fSSEnhLocMatSupport(NULL)
{
	
}

/* describe the parameters needed by the interface */
void GeoModelSS2DT::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  GeoModelSST::DefineParameters(list);
  
  /* 2D option must be plain strain */
  ParameterT& constraint = list.GetParameter("constraint_2D");
  constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void GeoModelSS2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	GeoModelSST::TakeParameterList(list);
  
	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusElas2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2)),
	fModulusContinuum2D.Dimension(dSymMatrixT::NumValues(2)),
	fModulusContinuumPerfPlas2D.Dimension(dSymMatrixT::NumValues(2)),
	fTotalStrain3D.Dimension(3);
}

/* returns elastic strain (3D) */
const dSymMatrixT& GeoModelSS2DT::ElasticStrain(const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return GeoModelSST::ElasticStrain(fTotalStrain3D, element, ip);
}


/* moduli */
const dMatrixT& GeoModelSS2DT::c_ijkl(void)
{
	/*SSSolidMat already reduces to 2D, trying to reduce again
	creates indexing errors */
	return SSSolidMatT::c_ijkl();

	/* 3D -> 2D, old method */
	/* fModulus2D.Rank4ReduceFrom3D(GeoModelSST::c_ijkl());
	return fModulus2D; */
}

const dMatrixT& GeoModelSS2DT::ce_ijkl(void)
{
	/* 3D -> 2D */
	fModulusElas2D.Rank4ReduceFrom3D(GeoModelSST::ce_ijkl());
	return fModulusElas2D;
}

const dMatrixT& GeoModelSS2DT::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulusPerfPlas2D.Rank4ReduceFrom3D(GeoModelSST::c_perfplas_ijkl());
	return fModulusPerfPlas2D;
}

const dMatrixT& GeoModelSS2DT::con_ijkl(void)
{
	/* 3D -> 2D */
	fModulusContinuum2D.Rank4ReduceFrom3D(GeoModelSST::con_ijkl());
	return fModulusContinuum2D;
}

const dMatrixT& GeoModelSS2DT::con_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulusContinuumPerfPlas2D.Rank4ReduceFrom3D(GeoModelSST::con_perfplas_ijkl());
	return fModulusContinuumPerfPlas2D;
}


/* stress */
const dSymMatrixT& GeoModelSS2DT::s_ij(void)
{
	fStress2D.ReduceFrom3D(GeoModelSST::s_ij());
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double GeoModelSS2DT::StrainEnergyDensity(void)
{
	return GeoModelSST::StrainEnergyDensity();
}
