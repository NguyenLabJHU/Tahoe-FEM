/* $Id: SolidMatList2DT.cpp,v 1.3 2001-04-27 10:53:29 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#include "SolidMatList2DT.h"

#include "ContinuumElementT.h"
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
#include "QuadLogOgden2DT.h"

/* constructor */
SolidMatList2DT::SolidMatList2DT(int length, const ElasticT& element_group):
	SolidMatListT(length),
	fElementGroup(element_group)
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
				fArray[matnum] = new SSKStV2D(in, fElementGroup);
				break;

			case kFDKStV:
				fArray[matnum] = new FDKStV2D(in, fElementGroup);
				break;

			case kSSCubic:
				fArray[matnum] = new SSCubic2DT(in, fElementGroup);
				break;

			case kFDCubic:
				fArray[matnum] = new FDCubic2DT(in, fElementGroup);
				break;

			case kSimoIso:
				fArray[matnum] = new SimoIso2D(in, fElementGroup);
				break;

			case kQuadLog:
				fArray[matnum] = new QuadLog2D(in, fElementGroup);
				break;

			case kQuadLogOgden:
				fArray[matnum] = new QuadLogOgden2DT(in, fElementGroup);												
				break;

			case kJ2SSKStV:
				fArray[matnum] = new J2SSKStV2D(in, fElementGroup);
				fHasHistory = true;															
				break;

			case kJ2Simo:
				fArray[matnum] = new J2Simo2D(in, fElementGroup);
				fHasHistory = true;
				break;

			case kJ2QL:
				fArray[matnum] = new J2QL2DLinHardT(in, fElementGroup);
				fHasHistory = true;
				break;

			case kDPSSKStV:
				fArray[matnum] = new DPSSKStV2D(in, fElementGroup);
				fHasHistory = true;															
				break;

			case kLJTr2D:
				fArray[matnum] = new LJTr2D(in, fElementGroup);
				break;

			case kLJFCC111:
				fArray[matnum] = new LJFCC111(in, fElementGroup);
				break;

			case kFCCEAM:
			{
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

				fArray[matnum] = new EAMFCC2D(in, fElementGroup, plane_code);			
				break;
			}

			case kmodCauchyBornDC:
			{
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

				fArray[matnum] = new ModCB2DT(in, fElementGroup, true, plane_code);
				break;
			}

			case kVIB:
				fArray[matnum] = new VIB2D(in, fElementGroup);
				fHasLocalizers = true;  				
				break;

			case kIsoVIBSimo:
				fArray[matnum] = new IsoVIB2D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIBOgden:
				fArray[matnum] = new OgdenIsoVIB2D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIBSimoJ2:
				fArray[matnum] = new J2IsoVIB2DLinHardT(in, fElementGroup);
				fHasLocalizers = true;
				fHasHistory = true;
				break;

			case kABAQUS_BCJ:
#ifdef __F2C__
				fArray[matnum] = new ABAQUS_BCJ(in, fElementGroup);
				fHasHistory = true;
#else
				cout << "\n SolidMatList2DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_BCJ << endl;
				throw eBadInputValue;
#endif /* __F2C__ */
				break;

//TEMP
#if 0
			case kSWDC100:
				fArray[matnum] = new SWDiamond100(in, fElementGroup);
				break;

			case kSWDC110:
				fArray[matnum] = new SWDiamond110(in, fElementGroup);
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
		int LTfnum = pmat->ThermalLTfNumber();
		if (LTfnum > -1)
		{
			pmat->SetThermalLTfPtr(fElementGroup.GetLTfPtr(LTfnum));
			
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
