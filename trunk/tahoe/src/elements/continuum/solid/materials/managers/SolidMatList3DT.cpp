/* $Id: SolidMatList3DT.cpp,v 1.20 2002-08-01 23:13:05 rdorgan Exp $ */
/* created: paklein (02/14/1997) */

#include "SolidMatList3DT.h"

#include "SmallStrainT.h"
#include "FiniteStrainT.h"
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
#include "SSStandardT.h"
#include "FDStandardT.h"
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
#include "tevp3D.h"
#include "LocalJ2SSNonlinHard.h"
#include "GradJ2SSNonlinHard.h"

#include "ABAQUS_BCJ.h"
#include "ABAQUS_VUMAT_BCJ.h"

//#include "OgdenViscVIB3D.h"

/* constructors */

using namespace Tahoe;

SolidMatList3DT::SolidMatList3DT(int length, const ElasticT& element_group):
	StructuralMatListT(length),
	fElementGroup(element_group)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatList3DT::SolidMatList3DT: WARNING: environment has no RTTI. Some\n" 
	     <<   "    consistency checking is disabled" << endl;
	/* cast and hope for the best */
	fSmallStrain = (const SmallStrainT*) &fElementGroup;
	fFiniteStrain = (const FiniteStrainT*) &fElementGroup;
#else

	/* cast to small strain */
	fSmallStrain = dynamic_cast<const SmallStrainT*>(&fElementGroup);

	/* cast to small strain */
	fFiniteStrain = dynamic_cast<const FiniteStrainT*>(&fElementGroup);
	
	/* must have at least one */
	if (!fSmallStrain && !fFiniteStrain)
	{
		cout << "\n SolidMatList3DT::SolidMatList3DT: could not cast element group to\n" 
		     <<   "     SmallStrainT or FiniteStrainT" << endl;
		throw eGeneralFail;
	}
#endif
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
		if (matnum < 0  || matnum >= fLength) throw eBadInputValue;

		/* repeated material number */
		if (fArray[matnum] != NULL)
		{
			cout << "\n SolidMatList3DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw eBadInputValue;
		}
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kLJTr2D:
			case kLJFCC111:
			{
				cout << "\n SolidMatList3DT::ReadMaterialData: material is 2D only: " 
				     << matcode << endl;
				throw eBadInputValue;
			}
			case kSSKStV:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new SSKStV(in, *fSmallStrain);
				break;
			}
			case kFDKStV:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDKStV(in, *fFiniteStrain);
				break;
			}							
			case kSSCubic:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new SSCubicT(in, *fSmallStrain);
				break;
			}
			case kFDCubic:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new FDCubicT(in, *fFiniteStrain);
				break;
			}
			case kSimoIso:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new SimoIso3D(in, *fFiniteStrain);
				break;
			}
			case kQuadLog:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new QuadLog3D(in, *fFiniteStrain);
				break;
			}
			case kQuadLogOgden:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new QuadLogOgden3DT(in, *fFiniteStrain);												
				break;
			}
			case kJ2SSKStV:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new J2SSKStV(in, *fSmallStrain);
				fHasHistory = true;														
				break;
			}
			case kJ2Simo:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2Simo3D(in, *fFiniteStrain);
				fHasHistory = true;														
				break;
			}
			case kJ2QL:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2QLLinHardT(in, *fFiniteStrain);
				fHasHistory = true;														
				break;
			}
			case kDPSSKStV:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new DPSSKStV(in, *fSmallStrain);
				fHasHistory = true;															
				break;
			}
			case kFCCEAM:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new EAMFCC3DMatT(in, *fFiniteStrain);
				break;
			}
			case kmodCauchyBornDC:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ModCB3DT(in, *fFiniteStrain, true);
				break;
			}
			case kVIB:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new VIB3D(in, *fFiniteStrain);
				fHasLocalizers = true;
				break;
			}
			case kOgdenViscVIB:
			{
			  /*				// check 
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenViscVIB3D(in, *fFiniteStrain);
				fHasLocalizers = true;
				fHasHistory = true;
			  */
			  cout << "\n Viscoelastic VIB model not yet implemented in 3D\n";
			  throw eGeneralFail;
			  break;
			}
			case kIsoVIBSimo:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new IsoVIB3D(in, *fFiniteStrain);
				fHasLocalizers = true;
				break;
			}
			case kIsoVIBOgden:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new OgdenIsoVIB3D(in, *fFiniteStrain);
				fHasLocalizers = true;
				break;
			}	
			case kSSStandard:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
						
				fArray[matnum] = new SSStandardT(in, *fSmallStrain);
				fHasHistory = true;
				break;				
			}
			case kFDStandard:
			{
				/* check */
				if (!fFiniteStrain) Error_no_small_strain(cout, matcode);
						
				fArray[matnum] = new FDStandardT(in, *fFiniteStrain);
				fHasHistory = true;
				break;				
			}
			case kIsoVIBSimoJ2:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new J2IsoVIB3DLinHardT(in, *fFiniteStrain);
				fHasLocalizers = true;
				fHasHistory = true;
				break;
			}	
			case kThermoViscoPlastic:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);
				
				fArray[matnum] = new tevp3D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kHyperEVP:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new HyperEVP3D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kBCJHypo:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypo3D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgKE:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageKE3D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kBCJHypoIsoDmgYC:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new BCJHypoIsoDamageYC3D(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
                        case kFDXtalElast:
                        {
                                /* check */
                                if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

                                fArray[matnum] = new FDCrystalElast(in, *fFiniteStrain);
                                break;
                        }
			case kLocXtalPlast:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlast_C:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlast_C(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kGrdXtalPlast:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlast(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlastFp:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlastFp(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocXtalPlastFp_C:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new LocalCrystalPlastFp_C(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kGrdXtalPlastFp:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new GradCrystalPlastFp(in, *fFiniteStrain);
				fHasHistory = true;
				break;
			}
			case kLocJ2SSNlHard:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new LocalJ2SSNonlinHard(in, *fSmallStrain);
				fHasHistory = true;														
				break;
			}
			case kGrdJ2SSNlHard:
			{
				/* check */
				if (!fSmallStrain) Error_no_small_strain(cout, matcode);
	
				fArray[matnum] = new GradJ2SSNonlinHard(in, *fSmallStrain);
				fHasHistory = true;														
				break;
			}
			case kABAQUS_BCJ:
			{
#ifdef __F2C__			
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_BCJ(in, *fFiniteStrain);
				fHasHistory = true;
#else
				cout << "\n SolidMatList3DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_BCJ << endl;
				throw eBadInputValue;
#endif /* __F2C__ */
	
				break;
			}			
			case kABAQUS_VUMAT_BCJ:
			{
#ifdef __F2C__			
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new ABAQUS_VUMAT_BCJ(in, *fFiniteStrain);
				fHasHistory = true;
#else
				cout << "\n SolidMatList3DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_VUMAT_BCJ << endl;
				throw eBadInputValue;
#endif /* __F2C__ */
	
				break;
			}			
//TEMP
#if 0
			case kIsoVIB_X:
			{
				/* check */
				if (!fFiniteStrain) Error_no_finite_strain(cout, matcode);

				fArray[matnum] = new IsoVIB3D_X(in, *fFiniteStrain);
				fHasLocalizers = true;
				break;
			}			
#endif //TEMP

			default:
			
				cout << "\n SolidMatList3DT::ReadMaterialData: unknown material code: ";
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
	throw eBadInputValue;
}

void SolidMatList3DT::Error_no_finite_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList3DT: material " << matcode
		<< " requires a finite strain element" << endl;
	throw eBadInputValue;
}
