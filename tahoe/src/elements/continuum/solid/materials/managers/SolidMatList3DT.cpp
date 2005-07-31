/* $Id: SolidMatList3DT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (02/14/1997)                                          */

#include "SolidMatList3DT.h"

#include "ElasticT.h"
#include "fstreamT.h"

/* 3D material type codes */
#include "SSKStV.h"
#include "FDKStV.h"
#include "SSCubicT.h"
#include "FDCubicT.h"
#include "QuadLog3D.h"
#include "VIB3D.h"
#include "IsoVIB3D.h"
#include "IsoVIB3D_X.h"
#include "J2IsoVIB3DLinHardT.h"
#include "ModCB3DT.h"
#include "EAMFCC3DMatT.h"
#include "SimoIso3D.h"
#include "DPSSKStV.h"
#include "J2SSKStV.h"
#include "J2QLLinHardT.h"
#include "OgdenIsoVIB3D.h"

#include "ABAQUS_BCJ.h"

/* 3D material type codes */
const int kSSKStV      = 1;
const int kFDKStV      = 2;
const int kSSCubic	   = 3;			
const int kFDCubic	   = 4;			
const int kVIB	       = 5;			
const int kIsoVIB      = 6;			
const int kQuadLog	   = 7;			
const int kIsoVIB_X	   = 8;			
const int kIsoVIB_J2   = 9;
const int kmodCBDC     = 10;
const int kCBDC        = 11; // no internal DOF
const int kEAM_FCC     = 12;
const int kSimoIso3D   = 13;
const int kJ2SSKStV    = 14;
const int kDPSSKStV    = 15;
const int kOgdenIsoVIB = 16;			
const int kJ2QL        = 18;

const int kABAQUS_BCJ = 80;

const int kMaterialMin = 1;
const int kMaterialMax = 100;

/* constructors */
SolidMatList3DT::SolidMatList3DT(int length, const ElasticT& element_group):
	SolidMatListT(length),
	fElementGroup(element_group)
{

}

/* read material data from the input stream */
void SolidMatList3DT::ReadMaterialData(ifstreamT& in)
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
			cout << "\n SolidMatList3DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw eBadInputValue;
		}
		
		/* add to the list of materials */
		switch (matcode)
		{
			case kSSKStV:
				fArray[matnum] = new SSKStV(in, fElementGroup);
				break;

			case kFDKStV:
				fArray[matnum] = new FDKStV(in, fElementGroup);
				break;
							
			case kSSCubic:
				fArray[matnum] = new SSCubicT(in, fElementGroup);
				break;

			case kFDCubic:
				fArray[matnum] = new FDCubicT(in, fElementGroup);
				break;

			case kQuadLog:
				fArray[matnum] = new QuadLog3D(in, fElementGroup);
				break;

			case kSimoIso3D:
				fArray[matnum] = new SimoIso3D(in, fElementGroup);
				break;
				
			case kVIB:
				fArray[matnum] = new VIB3D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIB:
				fArray[matnum] = new IsoVIB3D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kOgdenIsoVIB:
				fArray[matnum] = new OgdenIsoVIB3D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIB_X:
				fArray[matnum] = new IsoVIB3D_X(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIB_J2:
				fArray[matnum] = new J2IsoVIB3DLinHardT(in, fElementGroup);
				fHasLocalizers = true;
				fHasHistory = true;
				break;

			case kmodCBDC:
				fArray[matnum] = new ModCB3DT(in, fElementGroup, true);
				break;

			case kCBDC:
				fArray[matnum] = new ModCB3DT(in, fElementGroup, true);
				break;

			case kEAM_FCC:
				fArray[matnum] = new EAMFCC3DMatT(in, fElementGroup);
				break;

			case kDPSSKStV:
				fArray[matnum] = new DPSSKStV(in, fElementGroup);
				fHasHistory = true;															
				break;

			case kJ2SSKStV:
				fArray[matnum] = new J2SSKStV(in, fElementGroup);
				fHasHistory = true;														
				break;

			case kJ2QL:
				fArray[matnum] = new J2QLLinHardT(in, fElementGroup);
				fHasHistory = true;														
				break;

			case kABAQUS_BCJ:
#ifdef __F2C__			
				fArray[matnum] = new ABAQUS_BCJ(in, fElementGroup);
				fHasHistory = true;
#else
				cout << "\n SolidMatList3DT::ReadMaterialData: model requires f2c support: "
				     << kABAQUS_BCJ << endl;
				throw eBadInputValue;
#endif /* __F2C__ */
				break;

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
		cout << "\n SolidMatList3DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
}
