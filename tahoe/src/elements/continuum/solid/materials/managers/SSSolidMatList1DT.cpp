/* $Id: SSSolidMatList1DT.cpp,v 1.1.4.1 2004-04-08 07:33:04 paklein Exp $ */
#include "SSSolidMatList1DT.h"
#include "SSMatSupportT.h"
#include "fstreamT.h"

/* 1D material types codes */
/* Add small strain linear elastic material here */
#include "SSHookean1D.h"

#ifdef __DEVELOPMENT__
#include "DevelopmentElementsConfig.h"
#endif

#ifdef GRAD_SMALL_STRAIN_DEV
#include "GradJ2SS1D.h"
#include "J2SSKStV1D.h"
#endif

using namespace Tahoe;

/* constructor */
SSSolidMatList1DT::SSSolidMatList1DT(int length, const SSMatSupportT& support):
	SolidMatListT(length, support),
	fSSMatSupport(&support),
	fGradSSMatSupport(NULL)
{
#pragma message("check spatial dimension of material support")
	SetName("small_strain_material_1D");

#ifdef GRAD_SMALL_STRAIN_DEV
	/* cast to gradient enhanced small strain support */
	fGradSSMatSupport = TB_DYNAMIC_CAST(const GradSSMatSupportT*, fSSMatSupport);
#endif
}

SSSolidMatList1DT::SSSolidMatList1DT(void):
	fSSMatSupport(NULL),
	fGradSSMatSupport(NULL)	
{
	SetName("small_strain_material_1D");
}

/* read material data from the input stream */
void SSSolidMatList1DT::ReadMaterialData(ifstreamT& in)
{
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
		{
			cout << "\n SSSolidMatList1DT::ReadMaterialData: repeated material number: ";
			cout << matnum + 1 << endl;
			throw ExceptionT::kBadInputValue;
		}
		
		/* add to the list of matxxerials */
		switch (matcode)
		{
			case kSSKStV:
			{
				fArray[matnum] = new SSHookean1D(in, *fSSMatSupport);
				break;
		  	}
			case kGradJ2SS:
			{
#ifdef GRAD_SMALL_STRAIN_DEV
				fArray[matnum] = new GradJ2SS1D(in, *fGradSSMatSupport);
				fHasHistory = true;
				break;
#else
				ExceptionT::BadInputValue("SSSolidMatList1DT::ReadMaterialData", "GRAD_SMALL_STRAIN_DEV not enabled: %d", matcode);
#endif
			}
			case kJ2SSKStV1D:
			{
#ifdef GRAD_SMALL_STRAIN_DEV
				/* check */
				if (!fSSMatSupport) Error_no_small_strain(cout, matcode);
			
				fArray[matnum] = new J2SSKStV1D(in, *fSSMatSupport);
				fHasHistory = true;															
				break;
#else
				ExceptionT::BadInputValue("SSSolidMatList1DT::ReadMaterialData", "GRAD_SMALL_STRAIN_DEV not enabled: %d", matcode);
#endif
			}
			default:
			{
				cout << "\n SSSolidMatList1DT::ReadMaterialData: unknown material code: ";
				cout << matcode << '\n' << endl;
				throw ExceptionT::kBadInputValue;
			}
		}

		/* safe cast since all structural */
		SolidMaterialT* pmat = (SolidMaterialT*) fArray[matnum];
		/* verify construction */
		if (!pmat) throw ExceptionT::kOutOfMemory;
		
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
		cout << "\n SSSolidMatList1DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
}
