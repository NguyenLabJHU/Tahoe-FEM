/* $Id: SSSolidMatList3DT.cpp,v 1.1.2.1 2004-01-21 19:10:18 paklein Exp $ */
#include "SSSolidMatList3DT.h"
#include "SSMatSupportT.h"

#include "fstreamT.h"
#include "SolidMaterialsConfig.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentMaterialsConfig.h"
#endif

#include "SSKStV.h"
#include "SSCubicT.h"

#ifdef VISCOELASTICITY
#include "SSLinearVE3D.h"
#endif

#ifdef VISCOELASTIC_MATERIALS_DEV
#include "SSSV_KStV3D.h"
#endif

#ifdef J2PLASTICITY_MATERIALS_DEV
#include "SSJ2LinHardT.h"
#endif

#ifdef PLASTICITY_J2_MATERIAL
#include "J2SSKStV.h"
#include "LocalJ2SSNonlinHard.h"
#include "GradJ2SSNonlinHard.h"
#endif

#ifdef ABAQUS_MATERIAL
#ifdef ABAQUS_BCJ_MATERIAL_DEV
#include "ABAQUS_SS_BCJ_ISO.h"
#endif
#endif

#ifdef PLASTICITY_DP_MATERIAL
#include "DPSSKStV.h"
#endif

#ifdef FOSSUM_MATERIAL_DEV
#include "FossumSSIsoT.h"
#endif

#ifdef PLASTICITY_MR_MATERIAL_DEV
#include "MRSSKStV.h"
#endif

using namespace Tahoe;

/* constructors */
SSSolidMatList3DT::SSSolidMatList3DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support)
{
	SetName("small_strain_material_3D");
}

SSSolidMatList3DT::SSSolidMatList3DT(void):
	fSSMatSupport(NULL)
{
	SetName("small_strain_material_3D");
}

/* read material data from the input stream */
void SSSolidMatList3DT::ReadMaterialData(ifstreamT& in)
{
	const char caller[] = "SSSolidMatList3DT::ReadMaterialData";

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
				fArray[matnum] = new SSKStV(in, *fSSMatSupport);
				break;
			}
			case kSSCubic:
			{
				fArray[matnum] = new SSCubicT(in, *fSSMatSupport);
				break;
			}
			case kJ2SSKStV:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new J2SSKStV(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kDPSSKStV:
			{
#ifdef PLASTICITY_DP_MATERIAL
				fArray[matnum] = new DPSSKStV(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_DP_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kMRSSKStV:
			{
#ifdef PLASTICITY_MR_MATERIAL_DEV
				fArray[matnum] = new MRSSKStV(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_MR_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kFossumSSIso:
			{
#ifdef FOSSUM_MATERIAL_DEV
				fArray[matnum] = new FossumSSIsoT(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "FOSSUM_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kSSLinearVE:
			{
#ifdef VISCOELASTICITY
				fArray[matnum] = new SSLinearVE3D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTICITY not enabled: %d", matcode);
#endif
			}
			case kLocJ2SSNlHard:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new LocalJ2SSNonlinHard(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kGrdJ2SSNlHard:
			{
#ifdef PLASTICITY_J2_MATERIAL
				fArray[matnum] = new GradJ2SSNonlinHard(in, *fSSMatSupport);
				fHasHistory = true;														
				break;
#else
				ExceptionT::BadInputValue(caller, "PLASTICITY_J2_MATERIAL not enabled: %d", matcode);
#endif
			}
			case kABAQUS_BCJ_ISO:
			{
#ifdef __F2C__
#if defined(ABAQUS_MATERIAL) && defined(ABAQUS_BCJ_MATERIAL_DEV)
	
				/* small strain interface */
				fArray[matnum] = new ABAQUS_SS_BCJ_ISO(in, *fSSMatSupport);
	
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "ABAQUS_MATERIAL or ABAQUS_BCJ_MATERIAL_DEV not enabled: %d", matcode);
#endif
#else
				ExceptionT::BadInputValue(caller, "model requires f2c support: %d", kABAQUS_BCJ_ISO);
#endif /* __F2C__ */	
			}			
			case kSSSVKStV:
			{
#if VISCOELASTIC_MATERIALS_DEV
				fArray[matnum] = new SSSV_KStV3D(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "VISCOELASTIC_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			case kSSJ2LinHard:
			{
#if J2PLASTICITY_MATERIALS_DEV
				fArray[matnum] = new SSJ2LinHardT(in, *fSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue(caller, "J2PLASTICITY_MATERIALS_DEV not enabled: %d", matcode);
#endif
			}
			default:
				ExceptionT::BadInputValue(caller, "unknown material code: %d", matcode);
		}

		/* safe cast since all structural */
		SolidMaterialT* pmat = (SolidMaterialT*) fArray[matnum];

		/* verify construction */
		if (!pmat) ExceptionT::OutOfMemory(caller);
		
		/* set thermal LTf pointer */
		int LTfnum = pmat->ThermalStrainSchedule();
		if (LTfnum > -1)
		{
			pmat->SetThermalSchedule(fSolidMatSupport->Schedule(LTfnum));
			
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
