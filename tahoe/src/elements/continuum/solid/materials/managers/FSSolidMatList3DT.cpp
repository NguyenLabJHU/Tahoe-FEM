/* $Id: FSSolidMatList3DT.cpp,v 1.1.4.10 2004-06-19 23:28:09 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "FSSolidMatList3DT.h"

#include "fstreamT.h"
#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "FDKStV.h"
#include "FDCubicT.h"
#include "QuadLog3D.h"
#include "SimoIso3D.h"
#include "QuadLogOgden3DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC3DMatT.h"
#include "FCC3D.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB3DT.h"
#endif

#ifdef VIB_MATERIAL
#include "VIB3D.h"
#include "IsoVIB3D.h"
#include "J2IsoVIB3DLinHardT.h"
#include "OgdenIsoVIB3D.h"
#endif

#ifdef VISCOELASTICITY
#include "RGSplitT.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "FDSV_KStV3D.h"
#include "RGSplit3D.h"
#endif

#ifdef ELASTIC_OGDEN_MATERIAL_DEV
#include "OgdenMaterialT.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast.h"
#include "LocalCrystalPlast_C.h"
#include "GradCrystalPlast.h"
#include "LocalCrystalPlastFp.h"
#include "LocalCrystalPlastFp_C.h"
#include "GradCrystalPlastFp.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP3D.h"
#include "BCJHypo3D.h"
#include "BCJHypoIsoDamageKE3D.h"
#include "BCJHypoIsoDamageYC3D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2Simo3D.h"
#include "J2QLLinHardT.h"
#endif

#ifdef ELASTICITY_CRYSTAL_MATERIAL
#include "FDCrystalElast.h"
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp3D.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_BCJ.h"
#include "ABAQUS_BCJ_ISO.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif
#endif

#ifdef SIERRA_MATERIAL
#include "SIERRA_HypoElasticT.h"
#ifdef __FOSSUM__
#include "SIERRA_Isotropic_Geomaterial.h"
#endif /* __FOSSUM__ */
#endif /* SIERRA_MATERIAL */

using namespace Tahoe;

/* constructors */
FSSolidMatList3DT::FSSolidMatList3DT(int length, const FSMatSupportT& support):
	SolidMatListT(length, support),
	fFSMatSupport(&support)
{
#pragma message("check spatial dimension of material support")

	SetName("large_strain_material_3D");
}

FSSolidMatList3DT::FSSolidMatList3DT(void):
	fFSMatSupport(NULL)
{
	SetName("large_strain_material_3D");
}

/* read material data from the input stream */
void FSSolidMatList3DT::ReadMaterialData(ifstreamT& in)
{
ExceptionT::GeneralFail("FSSolidMatList3DT::ReadMaterialData");
#if 0
	const char caller[] = "FSSolidMatList3DT::ReadMaterialData";

	int i, matnum;
	SolidT::TypeT matcode;
	try {

	/* read material data */
	for (i = 0; i < fLength; i++)
	{
		in >> matnum; matnum--;
		in >> matcode;
		
		/* checks */
		if (matnum < 0  || matnum >= fLength) throw ExceptionT::kBadInputValue;

		/* repeated material number */
		if (fArray[matnum] != NULL)
			ExceptionT::BadInputValue(caller, "repeated material number: %d", matnum + 1);
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kFDKStV:
			{
				fArray[matnum] = new FDKStV(in, *fFSMatSupport);
				break;
			}							
			case kFDCubic:
			{
				fArray[matnum] = new FDCubicT(in, *fFSMatSupport);
				break;
			}
			case kSimoIso:
			{
				fArray[matnum] = new SimoIso3D(in, *fFSMatSupport);
				break;
			}
			case kQuadLog:
			{
				fArray[matnum] = new QuadLog3D(in, *fFSMatSupport);
				break;
			}
			case kQuadLogOgden:
			{
				fArray[matnum] = new QuadLogOgden3DT(in, *fFSMatSupport);												
				break;
			}
			case kJ2Simo:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new J2Simo3D(in, *fFSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kJ2QL:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new J2QLLinHardT(in, *fFSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFCCEAM:
			{
#ifdef CAUCHY_BORN_MATERIAL
				fArray[matnum] = new EAMFCC3DMatT(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFCC:
			{
#ifdef CAUCHY_BORN_MATERIAL
				fArray[matnum] = new FCC3D(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kmodCauchyBornDC:
			{
#ifdef MODCBSW_MATERIAL
				fArray[matnum] = new ModCB3DT(in, *fFSMatSupport, true);
				break;
#else
				ExceptionT::BadInputValue(caller, "MODCBSW_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kVIB:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new VIB3D(in, *fFSMatSupport);
				fHasLocalizers = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kIsoVIBSimo:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new IsoVIB3D(in, *fFSMatSupport);
				fHasLocalizers = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kIsoVIBOgden:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new OgdenIsoVIB3D(in, *fFSMatSupport);
				fHasLocalizers = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}	
			case kIsoVIBSimoJ2:
			{
#ifdef VIB_MATERIAL
				fArray[matnum] = new J2IsoVIB3DLinHardT(in, *fFSMatSupport);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}	
			case kRGSplitVE:
			{
#ifdef VISCOELASTICITY
				fArray[matnum] = new RGSplitT(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTICITY not enabled: %d", matcode);
#endif
			}
			case kThermoViscoPlastic:
			{
#ifdef THERMO_VISCO_PLASTIC_MATERIAL
				fArray[matnum] = new tevp3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "THERMO_VISCO_PLASTIC_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kHyperEVP:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new HyperEVP3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kBCJHypo:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new BCJHypo3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kBCJHypoIsoDmgKE:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new BCJHypoIsoDamageKE3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kBCJHypoIsoDmgYC:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				fArray[matnum] = new BCJHypoIsoDamageYC3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFDXtalElast:
			{
#ifdef ELASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new FDCrystalElast(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "ELASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new LocalCrystalPlast(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast_C:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new LocalCrystalPlast_C(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kGrdXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new GradCrystalPlast(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlastFp:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new LocalCrystalPlastFp(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlastFp_C:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new LocalCrystalPlastFp_C(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kGrdXtalPlastFp:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				fArray[matnum] = new GradCrystalPlastFp(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kABAQUS_BCJ:
			{
#ifdef __F2C__
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
				fArray[matnum] = new ABAQUS_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ);
#endif /* __F2C__ */	
			}			
			case kABAQUS_BCJ_ISO:
			{
#ifdef __F2C__
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
	
				/* small vs large strain elements */
				if (fFSMatSupport)
					fArray[matnum] = new ABAQUS_BCJ_ISO(in, *fFSMatSupport);
				else if (fSSMatSupport)
					fArray[matnum] = new ABAQUS_SS_BCJ_ISO(in, *fSSMatSupport);
				else
					ExceptionT::GeneralFail(caller);
					
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ_ISO);
#endif /* __F2C__ */	
			}			
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__			
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_VUMAT_BCJ);
#endif /* __F2C__ */
			}			
			case kRGSplit:
			{
#ifdef VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RGSplit3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kFDSVKStV:
			{
#if VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDSV_KStV3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSSSVKStV:
			{
#if VISCOELASTIC_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV3D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kOgdenMat:
			{
#if ELASTIC_OGDEN_MATERIAL_DEV
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenMaterialT(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "ELASTIC_OGDEN_MATERIAL_DEV not enabled: %d", matcode);
#endif
			}
			case kSSJ2LinHard:
			{
#if J2PLASTICITY_MATERIALS_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSJ2LinHardT(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "J2PLASTICITY_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSIERRA_Hypoelastic:
			{
#if SIERRA_MATERIAL
				fArray[matnum] = new SIERRA_HypoElasticT(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIERRA_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSIERRA_Iso_Geomat:
			{
#if SIERRA_MATERIAL
#ifdef __FOSSUM__
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SIERRA_Isotropic_Geomaterial(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else /* __FOSSUM__ */
				ExceptionT::BadInputValue(caller, "requires module fossum");
#endif /* __FOSSUM__ */
#else /* SIERRA_MATERIAL */
				ExceptionT::BadInputValue(caller, "SIERRA_MATERIAL not enabled: %d", matcode);
#endif /* SIERRA_MATERIAL*/
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
		}

		/* safe cast since all structural */
		SolidMaterialT* pmat = (SolidMaterialT*) fArray[matnum];

		/* verify construction */
		if (!pmat) ExceptionT::OutOfMemory(caller);
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fSolidMatSupport->Schedule(LTfnum));
			
			/* set flag */
			fHasThermal = true;
		}				

		/* perform initialization */
		pmat->Initialize();			
	}  } /* end try */

	catch (ExceptionT::CodeT error)
	{
		ExceptionT::Throw(error, caller, "exception constructing material %d, index %d, code %d",
			i+1, matnum+1, matcode);
	}
#endif
}

/* information about subordinate parameter lists */
void FSSolidMatList3DT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	SolidMatListT::DefineSubs(sub_list);
	
	/* choice of 2D materials */
	sub_list.AddSub("fs_material_list_3D", ParameterListT::Once, true);
}

/* return the description of the given inline subordinate parameter list */
void FSSolidMatList3DT::DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		SubListT& sub_sub_list) const
{
	if (sub == "fs_material_list_3D")
	{
		order = ParameterListT::Choice;
	
		sub_sub_list.AddSub("large_strain_cubic");
		sub_sub_list.AddSub("large_strain_StVenant");
		sub_sub_list.AddSub("Simo_isotropic");
		sub_sub_list.AddSub("quad_log");
		sub_sub_list.AddSub("quad_log_Ogden");

#ifdef PLASTICITY_J2_MATERIAL
		sub_sub_list.AddSub("Simo_J2");
		sub_sub_list.AddSub("quad_log_J2");
#endif

#ifdef CAUCHY_BORN_MATERIAL
		sub_sub_list.AddSub("FCC_3D");
		sub_sub_list.AddSub("FCC_EAM");		
#endif

#ifdef MODCBSW_MATERIAL
		sub_sub_list.AddSub("Cauchy-Born_diamond");
#endif

#ifdef VIB_MATERIAL
		sub_sub_list.AddSub("VIB");
		sub_sub_list.AddSub("isotropic_VIB");
		sub_sub_list.AddSub("Ogden_isotropic_VIB");
#endif
	}
	else /* inherited */
		SolidMatListT::DefineInlineSub(sub, order, sub_sub_list);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* FSSolidMatList3DT::NewSub(const StringT& list_name) const
{
	/* try to construct material */
	FSSolidMatT* fs_solid_mat = NewFSSolidMat(list_name);
	if (fs_solid_mat)
		return fs_solid_mat;
	else /* inherited */
		return SolidMatListT::NewSub(list_name);
}

/* accept parameter list */
void FSSolidMatList3DT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	SolidMatListT::TakeParameterList(list);

	/* construct materials - NOTE: subs have been defined as a choice, but
	 * here we construct as many materials as are passed in */
	AutoArrayT<FSSolidMatT*> materials;
	const ArrayT<ParameterListT>& subs = list.Lists();
	for (int i = 0; i < subs.Length(); i++) {
		const ParameterListT& sub = subs[i];
		FSSolidMatT* mat = NewFSSolidMat(sub.Name());
		if (mat) {
			materials.Append(mat);
			mat->TakeParameterList(sub);

			/* set flags */
			if (mat->HasHistory()) fHasHistory = true;	
			if (mat->HasThermalStrain()) fHasThermal = true;
			if (mat->HasLocalization()) fHasLocalizers = true;
		}
	}

	/* transfer */
	Dimension(materials.Length());
	for (int i = 0; i < materials.Length(); i++) 
		fArray[i] = materials[i];	
}

/* construct the specified material or NULL if the request cannot be completed */
FSSolidMatT* FSSolidMatList3DT::NewFSSolidMat(const StringT& name) const
{
	FSSolidMatT* mat = NULL;

	if (name == "large_strain_cubic")
		mat = new FDCubicT;
	else if (name == "large_strain_StVenant")
		mat = new FDKStV;
	else if (name == "Simo_isotropic")
		mat = new SimoIso3D;
	else if (name == "quad_log")
		mat = new QuadLog3D;
	else if (name == "quad_log_Ogden")
		mat = new QuadLogOgden3DT;

#ifdef PLASTICITY_J2_MATERIAL
	else if (name == "Simo_J2")
		mat = new J2Simo3D;
	else if (name == "quad_log_J2")
		mat = new J2QLLinHardT;
#endif

#ifdef CAUCHY_BORN_MATERIAL
	else if (name == "FCC_3D")
		mat = new FCC3D;
	else if (name == "FCC_EAM")
		mat = new EAMFCC3DMatT;
#endif

#ifdef MODCBSW_MATERIAL
	else if (name == "Cauchy-Born_diamond")
		mat= new ModCB3DT;
#endif

#ifdef VIB_MATERIAL
	else if (name == "VIB")
		mat= new VIB3D;
	else if (name == "isotropic_VIB")
		mat= new IsoVIB3D;
	else if (name == "Ogden_isotropic_VIB")
		mat= new OgdenIsoVIB3D;
#endif

	/* set support */
	if (mat) mat->SetFSMatSupport(fFSMatSupport);

	return mat;

}
