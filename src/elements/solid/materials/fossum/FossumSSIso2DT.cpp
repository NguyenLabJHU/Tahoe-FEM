/* $Id: FossumSSIso2DT.cpp,v 1.15 2006-06-18 21:34:58 regueiro Exp $ */
#include "FossumSSIso2DT.h"

#include "SSEnhLocMatSupportT.h"

#include "ElementCardT.h"
#include "StringT.h"

#include "DevelopmentElementsConfig.h"

using namespace Tahoe;

/* constructor */
FossumSSIso2DT::FossumSSIso2DT(void):
	ParameterInterfaceT("Fossum_small_strain_2D"),
	fSSEnhLocMatSupport(NULL)
{
	
}

/* describe the parameters needed by the interface */
void FossumSSIso2DT::DefineParameters(ParameterListT& list) const
{
  /* inherited */
  FossumSSIsoT::DefineParameters(list);
  
  /* 2D option must be plain stress */
  ParameterT& constraint = list.GetParameter("constraint_2D");
  constraint.SetDefault(kPlaneStrain);
}

/* accept parameter list */
void FossumSSIso2DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FossumSSIsoT::TakeParameterList(list);
  
	/* dimension work space */
	fStress2D.Dimension(2);
	fModulus2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusElas2D.Dimension(dSymMatrixT::NumValues(2));
	fModulusPerfPlas2D.Dimension(dSymMatrixT::NumValues(2)),
	fModulusContinuum2D.Dimension(dSymMatrixT::NumValues(2)),
	fModulusContinuumPerfPlas2D.Dimension(dSymMatrixT::NumValues(2)),
	fTotalStrain3D.Dimension(3);
  
	/* cast to small strain embedded discontinuity material pointer */
	//fSSEnhLocMatSupport = TB_DYNAMIC_CAST(const SSEnhLocMatSupportT*, fSSMatSupport);
}

#if 0
/* a pointer to the ParameterInterfaceT of the given subordinate */

ParameterInterfaceT* FossumSSIso2DT::NewSub(const StringT& name) const
{
  if (name == "Fossum_small_strain_2D")
    return new FossumSSIso2DT();
  else
    {
      /* inherited */
      ParameterInterfaceT* params = SSIsotropicMatT::NewSub(name);
      if (params) 
	return params;
      else
	return HookeanMatT::NewSub(name);
    }
}
#endif

/* initialization */
#if 0
void FossumSSIso2DT::Initialize(void)
{
ExceptionT::GeneralFail("FossumSSIso2DT::Initialize", "out of date");

	/* inherited */
	HookeanMatT::Initialize();

}
#endif


/* returns elastic strain (3D) */
const dSymMatrixT& FossumSSIso2DT::ElasticStrain(const dSymMatrixT& totalstrain, 
		const ElementCardT& element, int ip)
{
	/* 2D -> 3D (plane strain) */
	fTotalStrain3D.ExpandFrom2D(totalstrain);

	/* inherited */
	return FossumSSIsoT::ElasticStrain(fTotalStrain3D, element, ip);
}

#if 0
/* print parameters */
void FossumSSIso2DT::Print(ostream& out) const
{
	/* inherited */
	FossumSSIsoT::Print(out);
	Material2DT::Print(out);
}

/* print name */
void FossumSSIso2DT::PrintName(ostream& out) const
{
	/* inherited */
	FossumSSIsoT::PrintName(out);
	out << " 2D\n";
}
#endif

/* moduli */
const dMatrixT& FossumSSIso2DT::c_ijkl(void)
{
	/* 3D -> 2D */
	fModulus2D.Rank4ReduceFrom3D(FossumSSIsoT::c_ijkl());
//	fModulus2D *= fThickness;
	return fModulus2D;
}

const dMatrixT& FossumSSIso2DT::ce_ijkl(void)
{
	/* 3D -> 2D */
	fModulusElas2D.Rank4ReduceFrom3D(FossumSSIsoT::ce_ijkl());
	return fModulusElas2D;
}

const dMatrixT& FossumSSIso2DT::c_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulusPerfPlas2D.Rank4ReduceFrom3D(FossumSSIsoT::c_perfplas_ijkl());
//	fModulusPerfPlas2D *= fThickness;
	return fModulusPerfPlas2D;
}

const dMatrixT& FossumSSIso2DT::con_ijkl(void)
{
	/* 3D -> 2D */
	fModulusContinuum2D.Rank4ReduceFrom3D(FossumSSIsoT::con_ijkl());
//	fModulusContinuum2D *= fThickness;
	return fModulusContinuum2D;
}

const dMatrixT& FossumSSIso2DT::con_perfplas_ijkl(void)
{
	/* 3D -> 2D */
	fModulusContinuumPerfPlas2D.Rank4ReduceFrom3D(FossumSSIsoT::con_perfplas_ijkl());
//	fModulusContinuumPerfPlas2D *= fThickness;
	return fModulusContinuumPerfPlas2D;
}


/* stress */
const dSymMatrixT& FossumSSIso2DT::s_ij(void)
{
#ifdef ENHANCED_STRAIN_LOC_DEV
	int ip = CurrIP();
	ElementCardT& element = CurrentElement();
	int elem = CurrElementNumber();
	int element_locflag = 0;
	if (element.IsAllocated()) 
	{
		element_locflag = fSSEnhLocMatSupport->ElementLocflag(elem);
	}
	if ( element_locflag == 2 )
	{
		fStress2D = fSSEnhLocMatSupport->ElementStress(elem,ip);
	}
	else
	{
		/* 3D -> 2D */
		fStress2D.ReduceFrom3D(FossumSSIsoT::s_ij());
	}
#else
	/* 3D -> 2D */
	fStress2D.ReduceFrom3D(FossumSSIsoT::s_ij());
#endif
	//	fStress2D *= fThickness; 
	return fStress2D;
}

/* returns the strain energy density for the specified strain */
double FossumSSIso2DT::StrainEnergyDensity(void)
{
	return FossumSSIsoT::StrainEnergyDensity();
}
