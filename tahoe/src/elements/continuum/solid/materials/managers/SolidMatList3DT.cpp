/* $Id: SolidMatList3DT.cpp,v 1.3 2001-04-27 10:53:29 paklein Exp $ */
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
#include "QuadLogOgden3DT.h"

#include "ABAQUS_BCJ.h"

/* constructors */
SolidMatList3DT::SolidMatList3DT(int length, const ElasticT& element_group):
	SolidMatListT(length),
	fElementGroup(element_group)
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
			case kJ2Simo:
			case kLJTr2D:
			case kLJFCC111:
				cout << "\n SolidMatList3DT::ReadMaterialData: model " << matcode
				     << " is not implemented in 3D" << endl;
				throw eBadInputValue;

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

			case kSimoIso:
				fArray[matnum] = new SimoIso3D(in, fElementGroup);
				break;

			case kQuadLog:
				fArray[matnum] = new QuadLog3D(in, fElementGroup);
				break;

			case kQuadLogOgden:
				fArray[matnum] = new QuadLogOgden3DT(in, fElementGroup);												
				break;

			case kJ2SSKStV:
				fArray[matnum] = new J2SSKStV(in, fElementGroup);
				fHasHistory = true;														
				break;

			case kJ2QL:
				fArray[matnum] = new J2QLLinHardT(in, fElementGroup);
				fHasHistory = true;														
				break;

			case kDPSSKStV:
				fArray[matnum] = new DPSSKStV(in, fElementGroup);
				fHasHistory = true;															
				break;

			case kFCCEAM:
				fArray[matnum] = new EAMFCC3DMatT(in, fElementGroup);
				break;

			case kmodCauchyBornDC:
				fArray[matnum] = new ModCB3DT(in, fElementGroup, true);
				break;

			case kVIB:
				fArray[matnum] = new VIB3D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIBSimo:
				fArray[matnum] = new IsoVIB3D(in, fElementGroup);
				fHasLocalizers = true;
				break;

			case kIsoVIBOgden:
				fArray[matnum] = new OgdenIsoVIB3D(in, fElementGroup);
				fHasLocalizers = true;
				break;				

			case kIsoVIBSimoJ2:
				fArray[matnum] = new J2IsoVIB3DLinHardT(in, fElementGroup);
				fHasLocalizers = true;
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

//TEMP
#if 0

			case kIsoVIB_X:
				fArray[matnum] = new IsoVIB3D_X(in, fElementGroup);
				fHasLocalizers = true;
				break;

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
