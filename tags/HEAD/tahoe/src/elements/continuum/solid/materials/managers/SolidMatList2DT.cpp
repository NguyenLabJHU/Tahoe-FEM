/* $Id: SolidMatList2DT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
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

/* 2D material types codes */
const int kSSKStV          = 1;			
const int kFDKStV          = 2;			
const int kSSCubic         = 3;			
const int kFDCubic         = 4;			
const int kLJTr2D          = 5;
const int kLJFCC111        = 6;
const int kSWDC100         = 7;	//improper CB material
const int kSWDC110         = 8;	//improper CB material
const int kEAMFCC2D100     = 9;
const int kEAMFCC2D110     = 10;
const int kEAMFCC2D111     = 11;
const int kVIB             = 12;
const int kModCB           = 13;
const int kSimoIso2D       = 14;
const int kJ2SSKStV2D      = 15;
const int kJ2Simo2D        = 16;
const int kQuadLog2D       = 17;
const int kJ2QL2D          = 18;
const int kIsoVIB          = 19;
const int kOgdenIsoVIB     = 20;
const int kJ2VIB2D         = 21; // plane strain VIB 3D with J2
const int kDPSSKStV2D      = 22;
const int kD2VIB           = 23; // plane stress VIB + gradient terms

const int kABAQUS_BCJ      = 80;

const int kMaterialMin = 1;
const int kMaterialMax = 100;

/* constructor */
SolidMatList2DT::SolidMatList2DT(int length, const ElasticT& element_group):
	SolidMatListT(length),
	fElementGroup(element_group)
{

}

/* read material data from the input stream */
void SolidMatList2DT::ReadMaterialData(ifstreamT& in)
{
	int i, matnum, matcode;
	try {

	/* read material data */
	for (i = 0; i < fLength; i++)
	{
		in >> matnum; matnum--;
		in >> matcode;
		
		/* checks */
		if (matnum < 0  || matnum >= fLength) throw eBadInputValue;
		if (matcode < kMaterialMin ||
		    matcode > kMaterialMax) throw eBadInputValue;
		
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

			case kLJTr2D:
				fArray[matnum] = new LJTr2D(in, fElementGroup);
				break;

			case kLJFCC111:
				fArray[matnum] = new LJFCC111(in, fElementGroup);
				break;

			case kSWDC100:
				fArray[matnum] = new SWDiamond100(in, fElementGroup);
				break;

			case kSWDC110:
				fArray[matnum] = new SWDiamond110(in, fElementGroup);
				break;

			case kEAMFCC2D100:
				fArray[matnum] = new EAMFCC2D(in, fElementGroup, EAMFCC2D::kFCC2Dnatural);			
				break;

			case kEAMFCC2D110:
				fArray[matnum] = new EAMFCC2D(in, fElementGroup, EAMFCC2D::kFCC2D110);			
				break;

			case kEAMFCC2D111:
				fArray[matnum] = new EAMFCC2D(in, fElementGroup, EAMFCC2D::kFCC2D111);			
				break;

			case kVIB:
				fArray[matnum] = new VIB2D(in, fElementGroup);
				fHasLocalizers = true;  				
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
			case kModCB:
				fArray[matnum] = new ModCB2DT(in, fElementGroup, true);
				break;

			case kSimoIso2D:
				fArray[matnum] = new SimoIso2D(in, fElementGroup);
				break;

			case kJ2SSKStV2D:
				fArray[matnum] = new J2SSKStV2D(in, fElementGroup);
				fHasHistory = true;															
				break;
				
			case kDPSSKStV2D:
				fArray[matnum] = new DPSSKStV2D(in, fElementGroup);
				fHasHistory = true;															
				break;

			case kJ2Simo2D:
				fArray[matnum] = new J2Simo2D(in, fElementGroup);
				fHasHistory = true;
				break;

			case kQuadLog2D:
				fArray[matnum] = new QuadLog2D(in, fElementGroup);
				break;

			case kJ2QL2D:
				fArray[matnum] = new J2QL2DLinHardT(in, fElementGroup);
				fHasHistory = true;
				break;

			case kJ2VIB2D:
				fArray[matnum] = new J2IsoVIB2DLinHardT(in, fElementGroup);
				fHasLocalizers = true;
				fHasHistory = true;
				break;

			case kIsoVIB:
				fArray[matnum] = new IsoVIB2D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kOgdenIsoVIB:
				fArray[matnum] = new OgdenIsoVIB2D(in, fElementGroup);
				fHasLocalizers = true;
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
