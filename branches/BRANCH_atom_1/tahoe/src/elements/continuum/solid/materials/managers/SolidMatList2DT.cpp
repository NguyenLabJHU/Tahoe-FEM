/* $Id: SolidMatList2DT.cpp,v 1.26 2002-11-14 17:06:22 paklein Exp $ */
/* created: paklein (02/14/1997) */
#include "SolidMatList2DT.h"
#include "fstreamT.h"

/* 2D material types codes */
#include "SSKStV2D.h"
#include "FDKStV2D.h"
#include "SSCubic2DT.h"
#include "FDCubic2DT.h"
#include "LJTr2D.h"
#include "LJFCC111.h"
#include "SWDiamond110.h"
#include "SWDiamond100.h"
#include "EAMFCC2D.h"
#include "VIB2D.h"
#include "ModCB2DT.h"
#include "SimoIso2D.h"
#include "J2SSKStV2D.h"
#include "J2Simo2D.h"
#include "QuadLog2D.h"
#include "FossumSSIso2DT.h"
#include "J2QL2DLinHardT.h"
#include "IsoVIB2D.h"
#include "J2IsoVIB2DLinHardT.h"
#include "VoterChenCu.h"
#include "DPSSKStV2D.h"
#include "D2VIB2D_a.h"
#include "OgdenIsoVIB2D.h"
#include "LocalJ2SSNonlinHard2D.h"
#include "GradJ2SSNonlinHard2D.h"
#include "ABAQUS_BCJ.h"
#include "ABAQUS_VUMAT_BCJ.h"
#include "QuadLogOgden2DT.h"
#include "RGVIB2D.h"
#include "RG_NeoHookean2D.h"
#include "SV_NeoHookean2D.h"
#include "SSSV_KStV2D.h"
#include "FDSV_KStV2D.h"
#include "tevp2D.h"
#include "povirk2D.h"

#include "HyperEVP2D.h"
#include "BCJHypo2D.h"
#include "BCJHypoIsoDamageKE2D.h"
#include "BCJHypoIsoDamageYC2D.h"
#include "LocalCrystalPlast2D.h"
#include "GradCrystalPlast2D.h"
#include "LocalCrystalPlastFp2D.h"
#include "GradCrystalPlastFp2D.h"

using namespace Tahoe;

/* constructor */
SolidMatList2DT::SolidMatList2DT(int length, const StructuralMatSupportT& support):
	StructuralMatListT(length, support)
{

}

/* read material data from the input stream */
void SolidMatList2DT::ReadMaterialData(ifstreamT& in)
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
			cout << "\n SolidMatList2DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw ExceptionT::kBadInputValue;
		}
		
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
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDKStV2D(in, *fFDMatSupport);
				break;
			}
			/* case kMSKStV: // ---- Multi-Scale Kirchhoff - St. Venant
			{
				// check 
				if (!fMultiScale) Error_no_multi_scale(cout, matcode);

				fArray[matnum] = new MSKStV2D(in, *fMultiScale);
				break;
			} */
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
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new FDCubic2DT(in, *fFDMatSupport);
				break;
			}
			case kSimoIso:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new SimoIso2D(in, *fFDMatSupport);
				break;
			}
			case kQuadLog:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new QuadLog2D(in, *fFDMatSupport);
				break;
			}
			case kQuadLogOgden:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new QuadLogOgden2DT(in, *fFDMatSupport);												
				break;
			}
			case kJ2SSKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new J2SSKStV2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
			}
			case kJ2Simo:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new J2Simo2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kJ2QL:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new J2QL2DLinHardT(in, *fFDMatSupport);
				fHasHistory = true;
				break;
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
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LJTr2D(in, *fFDMatSupport);
				break;
			}
			case kLJFCC111:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LJFCC111(in, *fFDMatSupport);
				break;
			}
			case kFCCEAM:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				int i_plane_code;
				in >> i_plane_code;
				EAMFCC2D::PlaneCodeT plane_code = (EAMFCC2D::PlaneCodeT) i_plane_code;
				if (plane_code != EAMFCC2D::kFCC001 &&
				    plane_code != EAMFCC2D::kFCC101 &&
				    plane_code != EAMFCC2D::kFCC111)
				{
					cout << "\n SolidMatList2DT::ReadMaterialData: unrecognized plane code: "
					     << i_plane_code << endl;
					throw ExceptionT::kBadInputValue;
				}

				fArray[matnum] = new EAMFCC2D(in, *fFDMatSupport, plane_code);			
				break;
			}

			case kmodCauchyBornDC:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				int i_plane_code;
				in >> i_plane_code;
				ModCB2DT::PlaneCodeT plane_code = (ModCB2DT::PlaneCodeT) i_plane_code;
				if (plane_code != ModCB2DT::kDC001 &&
				    plane_code != ModCB2DT::kDC101 &&
				    plane_code != ModCB2DT::kDC111)
				{
					cout << "\n SolidMatList2DT::ReadMaterialData: unrecognized plane code: "
					     << i_plane_code << endl;
					throw ExceptionT::kBadInputValue;
				}

				fArray[matnum] = new ModCB2DT(in, *fFDMatSupport, true, plane_code);
				break;
			}

			case kVIB:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new VIB2D(in, *fFDMatSupport);
				fHasLocalizers = true;  				
				break;
			}
			case kIsoVIBSimo:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode); 
			
				fArray[matnum] = new IsoVIB2D(in, *fFDMatSupport);
				fHasLocalizers = true;
				break;
			}
			case kIsoVIBOgden:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new OgdenIsoVIB2D(in, *fFDMatSupport);
				fHasLocalizers = true;
				break;
			}
			case kIsoVIBSimoJ2:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new J2IsoVIB2DLinHardT(in, *fFDMatSupport);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
			}
		        case kFossumSSIso:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new FossumSSIso2DT(in, *fSSMatSupport);
				fHasHistory = true;		       		      				
				break;
			}
	
			case kThermoViscoPlastic:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new tevp2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kPovirk2D:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new povirk2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kHyperEVP:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new HyperEVP2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kBCJHypo:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypo2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgKE:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageKE2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgYC:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageYC2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlast:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlast_C:
			{
				cout << "\n SolidMatList2DT::ReadMaterialData: model " << matcode
				     << " is not implemented in 2D" << endl;
				throw ExceptionT::kBadInputValue;
				break;
			}
			case kGrdXtalPlast:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlast2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlastFp:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlastFp2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kGrdXtalPlastFp:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlastFp2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kLocJ2SSNlHard:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new LocalJ2SSNonlinHard2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
			}
			case kGrdJ2SSNlHard:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new GradJ2SSNonlinHard2D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
			}
			case kABAQUS_BCJ:
			{
#ifdef __F2C__
				fArray[matnum] = new ABAQUS_BCJ(in, *fFDMatSupport);
				fHasHistory = true;
#else
				cout << "\n SolidMatList2DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_BCJ << endl;
				throw ExceptionT::kBadInputValue;
#endif /* __F2C__ */
				break;
			}
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__
				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFDMatSupport);
				fHasHistory = true;
#else
				cout << "\n SolidMatList2DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_VUMAT_BCJ << endl;
				throw ExceptionT::kBadInputValue;
#endif /* __F2C__ */
				break;
			}
			case kRGVIB:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RGVIB2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kRGNeoHookean:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new RG_NeoHookean2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kSVNeoHookean:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SV_NeoHookean2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kFDSVKStV:
			{
				/* check */
				if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDSV_KStV2D(in, *fFDMatSupport);
				fHasHistory = true;
				break;
			}
			case kSSSVKStV:
			{
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);

				fArray[matnum] = new SSSV_KStV2D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
			}
//TEMP
#if 0
			case kSWDC100:
				fArray[matnum] = new SWDiamond100(in, *fFDMatSupport);
				break;

			case kSWDC110:
				fArray[matnum] = new SWDiamond110(in, *fFDMatSupport);
				break;


			case kD2VIB:
			{
#ifdef __NO_RTTI__
				cout << "\n SolidMatList2DT::ReadMaterialData: need RTTI" << endl;
				throw ExceptionT::kBadInputValue;
#endif

				const D2MeshFreeFDElasticT* D2ElementGroup;
				D2ElementGroup = dynamic_cast<const D2MeshFreeFDElasticT*>(&fElementGroup);
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

			/*			case kOgdenViscVIB:
						{
						if (!fFDMatSupport) Error_no_finite_strain(cout, matcode);
						
						fArray[matnum] = new OgdenViscVIB2D(in, *fFDMatSupport);
						fHasHistory = true;
						break;
						}*/
			default:
			
				cout << "\n SolidMatList2DT::ReadMaterialData: unknown material code: ";
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
		cout << "\n SolidMatList2DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
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

/* errror messages */

void SolidMatList2DT::Error_no_small_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList2DT: material " << matcode
		<< " requires a small strain element" << endl;
	throw ExceptionT::kBadInputValue;
}

void SolidMatList2DT::Error_no_finite_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList2DT: material " << matcode
		<< " requires a finite strain element" << endl;
	throw ExceptionT::kBadInputValue;
}

void SolidMatList2DT::Error_no_multi_scale(ostream& out, int matcode) const
{
	out << "\n SolidMatList2DT: material " << matcode
		<< " requires a variational multi-scale element" << endl;
	throw ExceptionT::kBadInputValue;
}

