/* $Id: Traction_CardT.cpp,v 1.9 2004-10-20 21:20:32 paklein Exp $ */
/* created: paklein (05/29/1996) */
#include "Traction_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "dArray2DT.h"
#include "ScheduleT.h"
#include "ElementSupportT.h"

#include "ElementsConfig.h"
#ifdef SHAPE_FUNCTION_CLASSES
#include "DomainIntegrationT.h"
#endif

using namespace Tahoe;

Traction_CardT::CoordSystemT Traction_CardT::int2CoordSystemT(int i)
{
	if (i == kCartesian)
		return kCartesian;
	else if (i == kLocal)
		return kLocal;
	else
		ExceptionT::GeneralFail("Traction_CardT::int2CoordSystemT", 
			"could not translate %d", i);
		
	return kLocal;
}

/* constructor */
Traction_CardT::Traction_CardT(void):
	fElemNum(0),
	fFacetNum(0),
	fCoordSystem(kCartesian),
	fLTfPtr(NULL),
	fValues(LocalArrayT::kUnspecified),
	fTau(0.0)
{
#ifndef SHAPE_FUNCTION_CLASSES
	ExceptionT::GeneralFail("Traction_CardT::Traction_CardT", "SHAPE_FUNCTION_CLASSES not enabled");
#endif
}	

/* modifiers */
void Traction_CardT::SetValues(const ElementSupportT& support, int elem, int facet,
	int nLTf, CoordSystemT coord_sys, const iArrayT& locnodenums,
	const dArray2DT& valuesT)
{	
	fValues.Dimension(valuesT.MajorDim(), valuesT.MinorDim());

	/* set */
	fElemNum  = elem;
	fFacetNum = facet;
	fCoordSystem = coord_sys;
	fLocNodeNums = locnodenums;
	fValues.FromTranspose(valuesT);

	/* resolve the pointer to the LTf */
	fLTfPtr = support.Schedule(nLTf);
}	

/* return the traction value: (ndof x nnd) */
void Traction_CardT::CurrentValue(LocalArrayT& traction) const
{
	if (fabs(fTau) < kSmall)
		traction.SetToScaled(fLTfPtr->Value(), fValues);
	else
		traction.SetToScaled(fTau*fLTfPtr->Rate() + fLTfPtr->Value(), fValues);
}

/* write the standard header */
void Traction_CardT::WriteHeader(ostream& out, int ndof) const
{
	double* junk = NULL;
	int d_width = OutputWidth(out, junk);
	out << setw(kIntWidth) << "elem";
	out << setw(kIntWidth) << "facet";
	out << setw(kIntWidth) << "LTf";
	out << setw(kIntWidth) << "coord";
	for (int i = 0; i < ndof; i++)
		out << setw(d_width - 2) << "t[" << i+1 << "]";
	out << '\n';			
}

namespace Tahoe {

/* input operator for codes */
istream& operator>>(istream& in, Traction_CardT::CoordSystemT& code)
{
	int i_code;
	in >> i_code;

	/* resolve code */
	switch (i_code)
	{
		case Traction_CardT::kCartesian:
			code = Traction_CardT::kCartesian;
			break;
		case Traction_CardT::kLocal:
			code = Traction_CardT::kLocal;
			break;
		default:
			cout << "\n operator>>Traction_CardT::CoordSystemT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
	}

	return in;
}

}
