/* $Id: SolidMatList3DT.cpp,v 1.28 2003-01-30 00:43:39 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "SolidMatList3DT.h"
#include "SolidMaterialsConfig.h"
#include "fstreamT.h"

/* 3D material type codes */
#include "SSKStV.h"
#include "FDKStV.h"
#include "SSCubicT.h"
#include "FDCubicT.h"
#include "QuadLog3D.h"
#include "SimoIso3D.h"
#include "DPSSKStV.h"
#include "QuadLogOgden3DT.h"
#include "tevp3D.h"
#include "LocalJ2SSNonlinHard.h"
#include "GradJ2SSNonlinHard.h"
#include "ABAQUS_BCJ.h"
#include "ABAQUS_VUMAT_BCJ.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC3DMatT.h"
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

#ifdef SIMO_HOLZAPFEL_MATERIAL
#include "SV_NeoHookean3D.h"
#include "SSSV_KStV3D.h"
#include "FDSV_KStV3D.h"
#endif

#ifdef REESE_GOVINDJEE_MATERIAL
#include "RG_NeoHookean3D.h"
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
#include "J2SSKStV.h"
#include "J2QLLinHardT.h"
#endif

#ifdef ELASTICITY_CRYSTAL_MATERIAL
#include "FDCrystalElast.h"
#endif

#ifdef FOSSUM_MATERIAL
#include "FossumSSIsoT.h"
#endif

using namespace Tahoe;

/* constructors */
SolidMatList3DT::SolidMatList3DT(int length, const SolidMatSupportT& support):
	SolidMatListT(length, support)
{

}

/* read material data from the input stream */
void SolidMatList3DT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "SolidMatList3DT::ReadMaterialData";

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
			case kLJTr2D:
			case kLJFCC111:
			{
				ExceptionT::BadInputValue(caller, "material is 2D only: %d", matcode);
			}
			case kSSKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new SSKStV(in, *fSSMatSupport);
				break;
			}
			case kFDKStV:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDKStV(in, *fFSMatSupport);
				break;
			}							
			case kSSCubic:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new SSCubicT(in, *fSSMatSupport);
				break;
			}
			case kFDCubic:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDCubicT(in, *fFSMatSupport);
				break;
			}
			case kSimoIso:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SimoIso3D(in, *fFSMatSupport);
				break;
			}
			case kQuadLog:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new QuadLog3D(in, *fFSMatSupport);
				break;
			}
			case kQuadLogOgden:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new QuadLogOgden3DT(in, *fFSMatSupport);												
				break;
			}
			case kJ2SSKStV:
			{
#ifdef PLASTICITY_J2_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new J2SSKStV(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kJ2Simo:
			{
#ifdef PLASTICITY_J2_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2QLLinHardT(in, *fFSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kDPSSKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new DPSSKStV(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
			}
			case kFCCEAM:
			{
#ifdef CAUCHY_BORN_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new EAMFCC3DMatT(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kmodCauchyBornDC:
			{
#ifdef MODCBSW_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ModCB3DT(in, *fFSMatSupport, true);
				break;
#else
				ExceptionT::BadInputValue(caller, "MODCBSW_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kVIB:
			{
#ifdef VIB_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2IsoVIB3DLinHardT(in, *fFSMatSupport);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VIB_MATERIAL not enabled: %d", matcode);
#endif
			}	
			case kFossumSSIso:
			{
#ifdef FOSSUM_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new FossumSSIsoT(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "FOSSUM_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kThermoViscoPlastic:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
				
				fArray[matnum] = new tevp3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
			}
			case kHyperEVP:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDCrystalElast(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "ELASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

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
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlastFp(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocJ2SSNlHard:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new LocalJ2SSNonlinHard(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
			}
			case kGrdJ2SSNlHard:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new GradJ2SSNonlinHard(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
			}
			case kABAQUS_BCJ:
			{
#ifdef __F2C__			
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ);
#endif /* __F2C__ */	
			}			
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__			
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_VUMAT_BCJ);
#endif /* __F2C__ */
			}			
			case kRGNeoHookean:
			{
#ifdef REESE_GOVINDJEE_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RG_NeoHookean3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "REESE_GOVINDJEE_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSVNeoHookean:
			{
#if SIMO_HOLZAPFEL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SV_NeoHookean3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIMO_HOLZAPFEL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFDSVKStV:
			{
#if SIMO_HOLZAPFEL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDSV_KStV3D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIMO_HOLZAPFEL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSSSVKStV:
			{
#if SIMO_HOLZAPFEL_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV3D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIMO_HOLZAPFEL_MATERIAL not enabled: %d", matcode);
#endif
			}
//TEMP
#if 0
			case kIsoVIB_X:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new IsoVIB3D_X(in, *fFSMatSupport);
				fHasLocalizers = true;
				break;
			}			
#endif
//TEMP
			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
		}

		/* safe cast since all structural */
		SolidMaterialT* pmat = (SolidMaterialT*) fArray[matnum];

		/* verify construction */
		if (!pmat) throw ExceptionT::kOutOfMemory;
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fSolidMatSupport.Schedule(LTfnum));
			
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

/* error messages */
void SolidMatList3DT::Error_no_small_strain(ostream& out, int matcode) const
{
#pragma unused(out)
#pragma unused(matcode)
	ExceptionT::BadInputValue("SolidMatList3DT::Error_no_small_strain", 
		"material %d requires a small strain element", matcode);
}

void SolidMatList3DT::Error_no_finite_strain(ostream& out, int matcode) const
{
#pragma unused(out)
#pragma unused(matcode)
	ExceptionT::BadInputValue("SolidMatList3DT::Error_no_small_strain", 
		"material %d requires a finite strain element", matcode);
}
