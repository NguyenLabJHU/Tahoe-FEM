/* $Id: Traction_CardT.cpp,v 1.7 2004-06-17 07:14:05 paklein Exp $ */
/* created: paklein (05/29/1996) */
#include "Traction_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"

#include "ifstreamT.h"
#include "dArray2DT.h"
#include "ScheduleT.h"
#include "ElementSupportT.h"

#include "ElementsConfig.h"
#ifdef SHAPE_FUNCTION_CLASSES
#include "DomainIntegrationT.h"
#endif

using namespace Tahoe;

/* constructor */
Traction_CardT::Traction_CardT(void):
	fElemNum(0),
	fFacetNum(0),
	fCoordSystem(kCartesian),
	fLTfPtr(NULL),
	fValues(LocalArrayT::kUnspecified)
{
#ifndef SHAPE_FUNCTION_CLASSES
	ExceptionT::GeneralFail("Traction_CardT::Traction_CardT", "SHAPE_FUNCTION_CLASSES not enabled");
#endif
}	

/* modifiers */
#ifdef SHAPE_FUNCTION_CLASSES
void Traction_CardT::EchoValues(const ElementSupportT& support, const DomainIntegrationT& domain,
	int elem, int ndof, ifstreamT& in, ostream& out)
{
	/* parameters */
	int facet;
	int nLTf;
	CoordSystemT coord_sys;
	dArray2DT valuesT;

	/* read  parameters */
	in >> facet >> nLTf >> coord_sys;

	/* correct offset */
	facet--;
	nLTf--;
	
	/* configure */
	domain.NodesOnFacet(facet, fLocNodeNums);
	valuesT.Dimension(fLocNodeNums.Length(), ndof);

	/* read tractions */
	in >> valuesT;
	
	/* set and echo */
	EchoValues(support, elem, facet, nLTf, coord_sys, fLocNodeNums, valuesT, out);
}
#else
void Traction_CardT::EchoValues(const ElementSupportT&, const DomainIntegrationT&, int, int, ifstreamT&, ostream&) {}
#endif	

void Traction_CardT::EchoValues(const ElementSupportT& support, int elem, int facet,
	int nLTf, CoordSystemT coord_sys, const iArrayT& locnodenums,
	const dArray2DT& valuesT, ostream& out)
{	
	fValues.Dimension(valuesT.MajorDim(), valuesT.MinorDim());

	/* set */
	fElemNum  = elem;
	fFacetNum = facet;
	fCoordSystem = coord_sys;
	fLocNodeNums = locnodenums;
	fValues.FromTranspose(valuesT);

	/* echo */
	out << setw(kIntWidth) << fElemNum + 1;
	out << setw(kIntWidth) << fFacetNum + 1;
	out << setw(kIntWidth) << nLTf + 1;
	out << setw(kIntWidth) << fCoordSystem;
	for (int i = 0; i < fLocNodeNums.Length(); i++)
	{
		/* tab */
		if (i > 0) out << setw(5*kIntWidth) << " ";
		valuesT.PrintRow(i, out);
	}

	/* resolve the pointer to the LTf */
	fLTfPtr = support.Schedule(nLTf);
}	

/* return the traction value: (ndof x nnd) */
void Traction_CardT::CurrentValue(LocalArrayT& traction) const
{
	traction.SetToScaled(fLTfPtr->Value(), fValues);
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
