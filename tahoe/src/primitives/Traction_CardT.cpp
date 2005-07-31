/* $Id: Traction_CardT.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/29/1996)                                          */

#include "Traction_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"

#include "fstreamT.h"
#include "dArray2DT.h"
#include "LoadTime.h"
#include "DomainIntegrationT.h"
#include "FEManagerT.h"

/* constructor */
Traction_CardT::Traction_CardT(void):
	fElemNum(0),
	fFacetNum(0),
	fCoordSystem(kCartesian),
	fLTfPtr(NULL),
	fValues(LocalArrayT::kUnspecified)
{

}	

/* modifiers */
void Traction_CardT::EchoValues(const FEManagerT& theBoss, const DomainIntegrationT& domain,
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
	valuesT.Allocate(fLocNodeNums.Length(), ndof);

	/* read tractions */
	in >> valuesT;
	
	/* set and echo */
	EchoValues(theBoss, elem, facet, nLTf, coord_sys, fLocNodeNums, valuesT, out);
}	

void Traction_CardT::EchoValues(const FEManagerT& theBoss, int elem, int facet,
	int nLTf, CoordSystemT coord_sys, const iArrayT& locnodenums,
	const dArray2DT& valuesT, ostream& out)
{	
	fValues.Allocate(valuesT.MajorDim(), valuesT.MinorDim());

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
	fLTfPtr = theBoss.GetLTfPtr(nLTf);
}	

/* return the traction value: (ndof x nnd) */
void Traction_CardT::CurrentValue(LocalArrayT& traction) const
{
	traction.SetToScaled(fLTfPtr->LoadFactor(), fValues);
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
			throw eBadInputValue;	
	}

	return in;
}
