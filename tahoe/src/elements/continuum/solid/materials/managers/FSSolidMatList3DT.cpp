/* $Id: FSSolidMatList3DT.cpp,v 1.1.2.1 2004-01-21 19:10:18 paklein Exp $ */
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
}
