#include "SolidMatList1DT.h"

#include "MultiScaleT.h"
#include "SmallStrainT.h"
#include "FiniteStrainT.h"

#include "fstreamT.h"

/* 1D material types codes */
/* Add small strain linear elastic material here */
#include "SSHookean1D.h"

/* constructor */
SolidMatList1DT::SolidMatList1DT(int length, const ElasticT& element_group):
	StructuralMatListT(length),
	fElementGroup(element_group)
{
#ifdef __NO_RTTI__
	cout << "\n SolidMatList1DT::SolidMatList1DT: WARNING: environment has no RTTI. Some\n" 
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
		cout << "\n SolidMatList1DT::SolidMatList1DT: could not cast element group to\n" 
		     <<   "     either SmallStrainT, FiniteStrainT, or MultiScaleT" << endl;
		throw eGeneralFail;
	}
#endif
}

/* read material data from the input stream */
void SolidMatList1DT::ReadMaterialData(ifstreamT& in)
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
			cout << "\n SolidMatList1DT::ReadMaterialData: repeated material number: ";
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
			cout << "Got to case kSSHookean1D!" << endl;
                        fArray[matnum] = new SSHookean1D(in, *fSmallStrain);
			cout << "Constructed new SSHookean1D!" << endl;
                        break;
		  }
			default:
			
				cout << "\n SolidMatList1DT::ReadMaterialData: unknown material code: ";
				cout << matcode << '\n' << endl;
				throw eBadInputValue;
		}

		/* safe cast since all structural */
		StructuralMaterialT* pmat = (StructuralMaterialT*) fArray[matnum];
		cout << "Past StructuralMaterialT!" << endl;
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
		cout << "\n SolidMatList1DT::ReadMaterialData: exception constructing material " << i+1
		     << '\n' << "     index " << matnum+1 << ", code " << matcode << endl;
		throw error;
	}
}


/* error messages */

void SolidMatList1DT::Error_no_small_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList1DT: material " << matcode
		<< " requires a small strain element" << endl;
	throw eBadInputValue;
}

void SolidMatList1DT::Error_no_finite_strain(ostream& out, int matcode) const
{
	out << "\n SolidMatList1DT: material " << matcode
		<< " requires a finite strain element" << endl;
	throw eBadInputValue;
}

void SolidMatList1DT::Error_no_multi_scale(ostream& out, int matcode) const
{
	out << "\n SolidMatList1DT: material " << matcode
		<< " requires a variational multi-scale element" << endl;
	throw eBadInputValue;
}

