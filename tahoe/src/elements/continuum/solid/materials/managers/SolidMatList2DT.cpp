/* $Id: SolidMatList2DT.cpp,v 1.30 2003-01-30 00:43:39 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "SolidMatList2DT.h"
#include "SolidMaterialsConfig.h"
#include "fstreamT.h"

/* 2D material types codes */
#include "SSKStV2D.h"
#include "FDKStV2D.h"
#include "SSCubic2DT.h"
#include "FDCubic2DT.h"
#include "SimoIso2D.h"
#include "QuadLog2D.h"
#include "QuadLogOgden2DT.h"

#ifdef CAUCHY_BORN_MATERIAL
#include "EAMFCC2D.h"
#include "LJTr2D.h"
#include "LJFCC111.h"
#endif

#ifdef MODCBSW_MATERIAL
#include "ModCB2DT.h"
#include "SWDiamond110.h"
#include "SWDiamond100.h"
#endif

#ifdef VIB_MATERIAL
#include "IsoVIB2D.h"
#include "J2IsoVIB2DLinHardT.h"
#include "VIB2D.h"
#include "D2VIB2D_a.h"
#include "OgdenIsoVIB2D.h"
#endif

#ifdef PLASTICITY_MACRO_MATERIAL
#include "HyperEVP2D.h"
#include "BCJHypo2D.h"
#include "BCJHypoIsoDamageKE2D.h"
#include "BCJHypoIsoDamageYC2D.h"
#endif

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "LocalCrystalPlast2D.h"
#include "GradCrystalPlast2D.h"
#include "LocalCrystalPlastFp2D.h"
#include "GradCrystalPlastFp2D.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2SSKStV2D.h"
#include "J2Simo2D.h"
#include "J2QL2DLinHardT.h"
#include "LocalJ2SSNonlinHard2D.h"
#include "GradJ2SSNonlinHard2D.h"
#endif

#ifdef PLASTICITY_DP_MATERIAL
#include "DPSSKStV2D.h"
#endif

#ifdef REESE_GOVINDJEE_MATERIAL
#include "RGVIB2D.h"
#include "RG_NeoHookean2D.h"
#endif

#ifdef SIMO_HOLZAPFEL_MATERIAL
#include "SV_NeoHookean2D.h"
#include "SSSV_KStV2D.h"
#include "FDSV_KStV2D.h"
#endif

#ifdef ABAQUS_MATERIAL
#include "ABAQUS_BCJ.h"
#include "ABAQUS_VUMAT_BCJ.h"
#endif

#ifdef THERMO_VISCO_PLASTIC_MATERIAL
#include "tevp2D.h"
#include "povirk2D.h"
#endif

#ifdef FOSSUM_MATERIAL
#include "FossumSSIso2DT.h"
#endif

using namespace Tahoe;

/* constructor */
SolidMatList2DT::SolidMatList2DT(int length, const SolidMatSupportT& support):
	SolidMatListT(length, support)
{

}

/* read material data from the input stream */
void SolidMatList2DT::ReadMaterialData(ifstreamT& in)
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
			case kSSKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new SSKStV2D(in, *fSSMatSupport);
				break;
			}
			case kFDKStV:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDKStV2D(in, *fFSMatSupport);
				break;
			}
			case kSSCubic:
			{ 
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new SSCubic2DT(in, *fSSMatSupport);
				break;
			}
			case kFDCubic:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new FDCubic2DT(in, *fFSMatSupport);
				break;
			}
			case kSimoIso:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new SimoIso2D(in, *fFSMatSupport);
				break;
			}
			case kQuadLog:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new QuadLog2D(in, *fFSMatSupport);
				break;
			}
			case kQuadLogOgden:
			{
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new QuadLogOgden2DT(in, *fFSMatSupport);												
				break;
			}
			case kJ2SSKStV:
			{
#ifdef PLASTICITY_J2_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new J2SSKStV2D(in, *fSSMatSupport);
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
			
				fArray[matnum] = new J2Simo2D(in, *fFSMatSupport);
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
			
				fArray[matnum] = new J2QL2DLinHardT(in, *fFSMatSupport);
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
			
				fArray[matnum] = new DPSSKStV2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
			}
			case kLJTr2D:
			{
#ifdef CAUCHY_BORN_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LJTr2D(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLJFCC111:
			{
#ifdef CAUCHY_BORN_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LJFCC111(in, *fFSMatSupport);
				break;
#else
				ExceptionT::BadInputValue(caller, "CAUCHY_BORN_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFCCEAM:
			{
#ifdef CAUCHY_BORN_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				int i_plane_code;
				in >> i_plane_code;
				EAMFCC2D::PlaneCodeT plane_code = (EAMFCC2D::PlaneCodeT) i_plane_code;
				if (plane_code != EAMFCC2D::kFCC001 &&
				    plane_code != EAMFCC2D::kFCC101 &&
				    plane_code != EAMFCC2D::kFCC111)
				{
					ExceptionT::BadInputValue(caller, "unrecognized plane code: %d", i_plane_code);
				}

				fArray[matnum] = new EAMFCC2D(in, *fFSMatSupport, plane_code);			
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

				int i_plane_code;
				in >> i_plane_code;
				ModCB2DT::PlaneCodeT plane_code = (ModCB2DT::PlaneCodeT) i_plane_code;
				if (plane_code != ModCB2DT::kDC001 &&
				    plane_code != ModCB2DT::kDC101 &&
				    plane_code != ModCB2DT::kDC111)
				{
					ExceptionT::BadInputValue(caller, "unrecognized plane code: %d", i_plane_code);
				}

				fArray[matnum] = new ModCB2DT(in, *fFSMatSupport, true, plane_code);
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
			
				fArray[matnum] = new VIB2D(in, *fFSMatSupport);
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
			
				fArray[matnum] = new IsoVIB2D(in, *fFSMatSupport);
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
			
				fArray[matnum] = new OgdenIsoVIB2D(in, *fFSMatSupport);
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
			
				fArray[matnum] = new J2IsoVIB2DLinHardT(in, *fFSMatSupport);
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
			
				fArray[matnum] = new FossumSSIso2DT(in, *fSSMatSupport);
				fHasHistory = true;		       		      				
				break;
#else
				ExceptionT::BadInputValue(caller, "FOSSUM_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kThermoViscoPlastic:
			{
#ifdef THERMO_VISCO_PLASTIC_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new tevp2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "THERMO_VISCO_PLASTIC_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kPovirk2D:
			{
#ifdef THERMO_VISCO_PLASTIC_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new povirk2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "THERMO_VISCO_PLASTIC_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kHyperEVP:
			{
#ifdef PLASTICITY_MACRO_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new HyperEVP2D(in, *fFSMatSupport);
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

				fArray[matnum] = new BCJHypo2D(in, *fFSMatSupport);
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

				fArray[matnum] = new BCJHypoIsoDamageKE2D(in, *fFSMatSupport);
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

				fArray[matnum] = new BCJHypoIsoDamageYC2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MACRO_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocXtalPlast_C:
			{
				ExceptionT::BadInputValue(caller, "model %d is not implemented in 2D", matcode);
			}
			case kGrdXtalPlast:
			{
#ifdef PLASTICITY_CRYSTAL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlast2D(in, *fFSMatSupport);
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

				fArray[matnum] = new LocalCrystalPlastFp2D(in, *fFSMatSupport);
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

				fArray[matnum] = new GradCrystalPlastFp2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_CRYSTAL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kLocJ2SSNlHard:
			{
#ifdef PLASTICITY_J2_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new LocalJ2SSNonlinHard2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kGrdJ2SSNlHard:
			{
#ifdef PLASTICITY_J2_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new GradJ2SSNonlinHard2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kABAQUS_BCJ:
			{
#ifdef __F2C__
#ifdef ABAQUS_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ);
#endif /* __F2C__ */
			}
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__
#ifdef ABAQUS_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ);
#endif /* __F2C__ */
			}
			case kRGVIB:
			{
#ifdef REESE_GOVINDJEE_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RGVIB2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "REESE_GOVINDJEE_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kRGNeoHookean:
			{
#ifdef REESE_GOVINDJEE_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RG_NeoHookean2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "REESE_GOVINDJEE_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSVNeoHookean:
			{
#ifdef SIMO_HOLZAPFEL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SV_NeoHookean2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIMO_HOLZAPFEL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFDSVKStV:
			{
#ifdef SIMO_HOLZAPFEL_MATERIAL
				/* check */
				if (!fFSMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDSV_KStV2D(in, *fFSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIMO_HOLZAPFEL_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSSSVKStV:
			{
#ifdef SIMO_HOLZAPFEL_MATERIAL
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "SIMO_HOLZAPFEL_MATERIAL not enabled: %d", matcode);
#endif
			}
//TEMP
#if 0
			case kSWDC100:
				fArray[matnum] = new SWDiamond100(in, *fFSMatSupport);
				break;

			case kSWDC110:
				fArray[matnum] = new SWDiamond110(in, *fFSMatSupport);
				break;


			case kD2VIB:
			{
#ifdef __NO_RTTI__
				cout << "\n SolidMatList2DT::ReadMaterialData: need RTTI" << endl;
				throw ExceptionT::kBadInputValue;
#endif
				const D2MeshFreeFSSolidT* D2ElementGroup;
				D2ElementGroup = dynamic_cast<const D2MeshFreeFSSolidT*>(&fElementGroup);
				if (!D2ElementGroup)
				{
					cout << "\n SolidMatList2DT::ReadMaterialData: material " << kD2VIB << " requires\n"
					     <<   "     a higher order gradient element" << endl;
					throw ExceptionT::kBadInputValue;
				}
			
				fArray[matnum] = new D2VIB2D_a(in, *D2ElementGroup);
				fHasLocalizers = true;  				
				break;
			}

#endif //TEMP

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

/* return true if the list contains plane stress models */
bool SolidMatList2DT::HasPlaneStress(void) const
{
	/* check materials */
	for (int i = 0; i < Length(); i++)
	{
		/* get pointer to Material2DT */
		const ContinuumMaterialT* cont_mat = fArray[i];
		const Material2DT* mat_2D = dynamic_cast<const Material2DT*>(cont_mat);
		if (!mat_2D) throw ExceptionT::kGeneralFail;

		/* test */
		if (mat_2D->ConstraintOption() == Material2DT::kPlaneStress) return true;
	}
	return false;
}

/* error messages */

void SolidMatList2DT::Error_no_small_strain(ostream& out, int matcode) const
{
#pragma unused(out)
#pragma unused(matcode)
	ExceptionT::BadInputValue("SolidMatList2DT::Error_no_small_strain", 
		"material %d requires a small strain element", matcode);
}

void SolidMatList2DT::Error_no_finite_strain(ostream& out, int matcode) const
{
#pragma unused(out)
#pragma unused(matcode)
	ExceptionT::BadInputValue("SolidMatList2DT::Error_no_small_strain", 
		"material %d requires a finite strain element", matcode);
}
