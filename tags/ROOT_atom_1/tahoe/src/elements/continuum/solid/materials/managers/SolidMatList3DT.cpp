/* $Id: SolidMatList3DT.cpp,v 1.24 2002-11-14 17:06:22 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "SolidMatList3DT.h"
#include "fstreamT.h"

/* 3D material type codes */
#include "SSKStV.h"
#include "FDKStV.h"
#include "SSCubicT.h"
#include "FDCubicT.h"
#include "QuadLog3D.h"
#include "VIB3D.h"
#include "IsoVIB3D.h"
#include "J2IsoVIB3DLinHardT.h"
#include "ModCB3DT.h"
#include "EAMFCC3DMatT.h"
#include "SimoIso3D.h"
#include "J2Simo3D.h"
#include "DPSSKStV.h"
#include "J2SSKStV.h"
#include "J2QLLinHardT.h"
#include "OgdenIsoVIB3D.h"
#include "QuadLogOgden3DT.h"
#include "FossumSSIsoT.h"
#include "HyperEVP3D.h"
#include "BCJHypo3D.h"
#include "BCJHypoIsoDamageKE3D.h"
#include "BCJHypoIsoDamageYC3D.h"
#include "FDCrystalElast.h"
#include "LocalCrystalPlast.h"
#include "LocalCrystalPlast_C.h"
#include "GradCrystalPlast.h"
#include "LocalCrystalPlastFp.h"
#include "LocalCrystalPlastFp_C.h"
#include "GradCrystalPlastFp.h"
#include "RG_NeoHookean3D.h"
#include "SV_NeoHookean3D.h"
#include "SSSV_KStV3D.h"
#include "FDSV_KStV3D.h"
#include "tevp3D.h"
#include "LocalJ2SSNonlinHard.h"
#include "GradJ2SSNonlinHard.h"

#include "ABAQUS_BCJ.h"
#include "ABAQUS_VUMAT_BCJ.h"

using namespace Tahoe;

/* constructors */
SolidMatList3DT::SolidMatList3DT(int length, const StructuralMatSupportT& support):
	StructuralMatListT(length, support)
{

}

/* read material data from the input stream */
void SolidMatList3DT::ReadMaterialData(ifstreamT& in)
{
	int i, matnum;
	MaterialT::SolidT matcode;
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
		{
			cout << "\n SolidMatList3DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw ExceptionT::kBadInputValue;
		}
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kLJTr2D:
			case kLJFCC111:
			{
				cout << "\n SolidMatList3DT::ReadMaterialData: material is 2D only: " 
				     << matcode << endl;
				throw ExceptionT::kBadInputValue;
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
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDKStV(in, *fFDMatSupport);
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
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDCubicT(in, *fFDMatSupport);
				break;
			}
			case kSimoIso:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SimoIso3D(in, *fFDMatSupport);
				break;
			}
			case kQuadLog:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new QuadLog3D(in, *fFDMatSupport);
				break;
			}
			case kQuadLogOgden:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new QuadLogOgden3DT(in, *fFDMatSupport);												
				break;
			}
			case kJ2SSKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new J2SSKStV(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
			}
			case kJ2Simo:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2Simo3D(in, *fFDMatSupport);
				fHasHistory = true;														
				break;
			}
			case kJ2QL:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2QLLinHardT(in, *fFDMatSupport);
				fHasHistory = true;														
				break;
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
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new EAMFCC3DMatT(in, *fFDMatSupport);
				break;
			}
			case kmodCauchyBornDC:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ModCB3DT(in, *fFDMatSupport, true);
				break;
			}
			case kVIB:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new VIB3D(in, *fFDMatSupport);
				fHasLocalizers = true;
				break;
			}
			/*			case kOgdenViscVIB:
						{
						if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
									       
						fArray[matnum] = new OgdenViscVIB3D(in, *fFDMatSupport);
						fHasLocalizers = true;
						fHasHistory = true;
						break;
						} */
			case kIsoVIBSimo:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new IsoVIB3D(in, *fFDMatSupport);
				fHasLocalizers = true;
				break;
			}
			case kIsoVIBOgden:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenIsoVIB3D(in, *fFDMatSupport);
				fHasLocalizers = true;
				break;
			}	
			case kIsoVIBSimoJ2:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2IsoVIB3DLinHardT(in, *fFDMatSupport);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
			}	
			case kFossumSSIso:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new FossumSSIsoT(in, *fSSMatSupport);
				fHasHistory = true;
				break;
			}
			case kThermoViscoPlastic:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
				
				fArray[matnum] = new tevp3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kHyperEVP:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new HyperEVP3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kBCJHypo:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypo3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgKE:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageKE3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgYC:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageYC3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kFDXtalElast:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDCrystalElast(in, *fFDMatSupport);
				break;
			}
			case kLocXtalPlast:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlast_C:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast_C(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kGrdXtalPlast:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlast(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlastFp:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlastFp(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlastFp_C:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlastFp_C(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kGrdXtalPlastFp:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlastFp(in, *fFDMatSupport);
				fHasHistory = true;
				break;
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
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_BCJ(in, *fFDMatSupport);
				fHasHistory = true;
#else
				cout << "\n SolidMatList3DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_BCJ << endl;
				throw ExceptionT::kBadInputValue;
#endif /* __F2C__ */
	
				break;
			}			
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__			
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFDMatSupport);
				fHasHistory = true;
#else
				cout << "\n SolidMatList3DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_VUMAT_BCJ << endl;
				throw ExceptionT::kBadInputValue;
#endif /* __F2C__ */
	
				break;
			}			
			case kRGNeoHookean:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RG_NeoHookean3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kSVNeoHookean:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SV_NeoHookean3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kFDSVKStV:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDSV_KStV3D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kSSSVKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV3D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
			}
//TEMP
#if 0
			case kIsoVIB_X:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new IsoVIB3D_X(in, *fFDMatSupport);
				fHasLocalizers = true;
				break;
			}			
#endif //TEMP

			default:
			
				cout << "\n SolidMatList3DT::ReadMaterialData: unknown material code: ";
				cout << matcode << '\n' << endl;
				throw ExceptionT::kBadInputValue;
		}

		/* safe cast since all structural */
		StructuralMaterialT* pmat = (StructuralMaterialT*) fArray[matnum];

		/* verify construction */
		if (!pmat) throw ExceptionT::kOutOfMemory;
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fStructuralMatSupport.Schedule(LTfnum));
			
			/* set flag */
			fHasThermal = true;
		}				

		/* perform initialization */
		pmat->Initialize();			
	}  } /* end try */

	catch (ExceptionT::CodeT error)
	{
		cout << "\n SolidMatList3DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
}

/* errror messages */
void SolidMatList3DT::Error_no_small_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList3DT: material " << matcode
		<< " requires a small strain element" << endl;
	throw ExceptionT::kBadInputValue;
}

void SolidMatList3DT::Error_no_finite_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList3DT: material " << matcode
		<< " requires a finite strain element" << endl;
	throw ExceptionT::kBadInputValue;
}
