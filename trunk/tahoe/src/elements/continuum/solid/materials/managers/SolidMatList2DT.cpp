/* $Id: SolidMatList2DT.cpp,v 1.18 2002-06-08 20:20:40 paklein Exp $ */
/* created: paklein (02/14/1997) */

#include "SolidMatList2DT.h"

#include "MultiScaleT.h"
#include "SmallStrainT.h"
#include "FiniteStrainT.h"
#include "D2MeshFreeFDElasticT.h"

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
#include "J2QL2DLinHardT.h"
#include "IsoVIB2D.h"
#include "J2IsoVIB2DLinHardT.h"
#include "VoterChenCu.h"
#include "DPSSKStV2D.h"
#include "D2VIB2D_a.h"
#include "OgdenIsoVIB2D.h"
#include "ABAQUS_BCJ.h"
#include "ABAQUS_VUMAT_BCJ.h"
#include "QuadLogOgden2DT.h"
#include "OgdenViscVIB2D.h"
#include "SKStVT2D.h"
#include "MaxwellT2D.h"
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

//#include "OgdenViscVIB2Dold.h"

/* constructor */
SolidMatList2DT::SolidMatList2DT(int length, const ElasticT& element_group):
	StructuralMatListT(length),
	fElementGroup(element_group)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatList2DT::SolidMatList2DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
	/* cast and hope for the best */
	fSmallStrain  = (const SmallStrainT*)  &fElementGroup;
	fFiniteStrain = (const FiniteStrainT*) &fElementGroup;
	fMultiScale   = (const MultiScaleT*)   &fElementGroup;
#else

	/* cast to small strain */
	fSmallStrain  = dynamic_cast<const SmallStrainT*>(&fElementGroup);

	/* cast to finite strain */
	fFiniteStrain = dynamic_cast<const FiniteStrainT*>(&fElementGroup);
	
	/* cast to multi scale */
	fMultiScale   = dynamic_cast<const MultiScaleT*>(&fElementGroup);
	
	/* must have at least one */
	if (!fSmallStrain && !fFiniteStrain && !fMultiScale)
	{
		cout << "\n SolidMatList2DT::SolidMatList2DT: could not cast element group to\n" 
		     <<   "     either SmallStrainT, FiniteStrainT, or MultiScaleT" << endl;
		throw eGeneralFail;
	}
#endif
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
		if (matnum < 0  || matnum >= fLength) throw eBadInputValue;
		
		/* repeated material number */
		if (fArray[matnum] != NULL)
		{
			cout << "\n SolidMatList2DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw eBadInputValue;
		}
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kSSKStV:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new SSKStV2D(in, *fSmallStrain);
				break;
			}
			case kFDKStV:
			{
				/* check */
				if (!fFiniteStrain && !fMultiScale) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDKStV2D(in, *fFiniteStrain);
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
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new SSCubic2DT(in, *fSmallStrain);
				break;
			}
			case kFDCubic:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new FDCubic2DT(in, *fFiniteStrain);
				break;
			}
			case kSimoIso:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new SimoIso2D(in, *fFiniteStrain);
				break;
			}
			case kQuadLog:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new QuadLog2D(in, *fFiniteStrain);
				break;
			}
			case kQuadLogOgden:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new QuadLogOgden2DT(in, *fFiniteStrain);												
				break;
			}
			case kJ2SSKStV:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new J2SSKStV2D(in, *fSmallStrain);
				fHasHistory = true;															
				break;
			}
			case kJ2Simo:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new J2Simo2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kJ2QL:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new J2QL2DLinHardT(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kDPSSKStV:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new DPSSKStV2D(in, *fSmallStrain);
				fHasHistory = true;															
				break;
			}
			case kLJTr2D:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LJTr2D(in, *fFiniteStrain);
				break;
			}
			case kLJFCC111:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LJFCC111(in, *fFiniteStrain);
				break;
			}
			case kFCCEAM:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				int i_plane_code;
				in >> i_plane_code;
				EAMFCC2D::PlaneCodeT plane_code = (EAMFCC2D::PlaneCodeT) i_plane_code;
				if (plane_code != EAMFCC2D::kFCC001 &&
				    plane_code != EAMFCC2D::kFCC101 &&
				    plane_code != EAMFCC2D::kFCC111)
				{
					cout << "\n SolidMatList2DT::ReadMaterialData: unrecognized plane code: "
					     << i_plane_code << endl;
					throw eBadInputValue;
				}

				fArray[matnum] = new EAMFCC2D(in, *fFiniteStrain, plane_code);			
				break;
			}

			case kmodCauchyBornDC:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				int i_plane_code;
				in >> i_plane_code;
				ModCB2DT::PlaneCodeT plane_code = (ModCB2DT::PlaneCodeT) i_plane_code;
				if (plane_code != ModCB2DT::kDC001 &&
				    plane_code != ModCB2DT::kDC101 &&
				    plane_code != ModCB2DT::kDC111)
				{
					cout << "\n SolidMatList2DT::ReadMaterialData: unrecognized plane code: "
					     << i_plane_code << endl;
					throw eBadInputValue;
				}

				fArray[matnum] = new ModCB2DT(in, *fFiniteStrain, true, plane_code);
				break;
			}

			case kVIB:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new VIB2D(in, *fFiniteStrain);
				fHasLocalizers = true;  				
				break;
			}
			case kIsoVIBSimo:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode); 
			
				fArray[matnum] = new IsoVIB2D(in, *fFiniteStrain);
				fHasLocalizers = true;
				break;
			}
			case kIsoVIBOgden:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new OgdenIsoVIB2D(in, *fFiniteStrain);
				fHasLocalizers = true;
				break;
			}
			case kIsoVIBSimoJ2:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new J2IsoVIB2DLinHardT(in, *fFiniteStrain);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
			}	
			case kThermoViscoPlastic:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new tevp2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
		        case kPovirk2D:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
			
				fArray[matnum] = new povirk2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kHyperEVP:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new HyperEVP2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kBCJHypo:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypo2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgKE:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageKE2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgYC:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageYC2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlast:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlast_C:
			{
				cout << "\n SolidMatList2DT::ReadMaterialData: model " << matcode
				     << " is not implemented in 2D" << endl;
				throw eBadInputValue;
				break;
			}
			case kGrdXtalPlast:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlast2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlastFp:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlastFp2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kGrdXtalPlastFp:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlastFp2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kABAQUS_BCJ:
			{
#ifdef __F2C__
				fArray[matnum] = new ABAQUS_BCJ(in, *fFiniteStrain);
				fHasHistory = true;
#else
				cout << "\n SolidMatList2DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_BCJ << endl;
				throw eBadInputValue;
#endif /* __F2C__ */
				break;
			}
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__
				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFiniteStrain);
				fHasHistory = true;
#else
				cout << "\n SolidMatList2DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_VUMAT_BCJ << endl;
				throw eBadInputValue;
#endif /* __F2C__ */
				break;
			}
			case kOgdenViscVIB:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenViscVIB2D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kSKStVT:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new SKStVT2D(in, *fSmallStrain);
				fHasHistory = true;
				break;
			}
			case kMaxwellT:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
				fArray[matnum] = new MaxwellT2D(in, *fSmallStrain);
				fHasHistory = true;
				break;
			}
//TEMP
#if 0
			case kSWDC100:
				fArray[matnum] = new SWDiamond100(in, *fFiniteStrain);
				break;

			case kSWDC110:
				fArray[matnum] = new SWDiamond110(in, *fFiniteStrain);
				break;


			case kD2VIB:
			{
#ifdef __NO_RTTI__
				cout << "\n SolidMatList2DT::ReadMaterialData: need RTTI" << endl;
				throw eBadInputValue;
#endif

				const D2MeshFreeFDElasticT* D2ElementGroup;
				D2ElementGroup = dynamic_cast<const D2MeshFreeFDElasticT*>(&fElementGroup);
				if (!D2ElementGroup)
				{
					cout << "\n SolidMatList2DT::ReadMaterialData: material " << kD2VIB << " requires\n"
					     <<   "     a higher order gradient element" << endl;
					throw eBadInputValue;
				}
			
				fArray[matnum] = new D2VIB2D_a(in, *D2ElementGroup);
				fHasLocalizers = true;  				
				break;
			}

#endif //TEMP

#if 0
			case kOgdenViscVIBold:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenViscVIB2Dold(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
#endif
			default:
			
				cout << "\n SolidMatList2DT::ReadMaterialData: unknown material code: ";
				cout << matcode << '\n' << endl;
				throw eBadInputValue;
		}

		/* safe cast since all structural */
		StructuralMaterialT* pmat = (StructuralMaterialT*) fArray[matnum];

		/* verify construction */
		if (!pmat) throw eOutOfMemory;
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fElementGroup.Schedule(LTfnum));
			
			/* set flag */
			fHasThermal = true;
		}
		
		/* perform initialization */
		pmat->Initialize();
	}  } /* end try */
	
	catch (int error)
	{
		cout << "\n SolidMatList2DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
}


/* errror messages */

void SolidMatList2DT::Error_no_small_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList2DT: material " << matcode
		<< " requires a small strain element" << endl;
	throw eBadInputValue;
}

void SolidMatList2DT::Error_no_finite_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList2DT: material " << matcode
		<< " requires a finite strain element" << endl;
	throw eBadInputValue;
}

void SolidMatList2DT::Error_no_multi_scale(ostream& out, int matcode) const
{
	out << "\n SolidMatList2DT: material " << matcode
		<< " requires a variational multi-scale element" << endl;
	throw eBadInputValue;
}

