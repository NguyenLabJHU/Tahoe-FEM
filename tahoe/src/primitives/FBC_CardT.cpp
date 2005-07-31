/* $Id: FBC_CardT.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (06/15/1996)                                          */

#include "FBC_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"

#include "fstreamT.h"
#include "NodeManagerT.h"
#include "LoadTime.h"

/* constructor */
FBC_CardT::FBC_CardT(void):
	fNode(-1),
	fDOF(-1),
	fLTf(-1),
	fValue(0.0),
	fLTfPtr(NULL)
{

}

/* modifiers */
void FBC_CardT::SetValues(const NodeManagerPrimitive& theBoss, ifstreamT& in)
{
	/* parameters */
	int    node;
	int    dof;
	int    nLTf;
	double value;

	/* read */
	in >> node >> dof >> nLTf >> value;
	
	/* correct offset */
	node--; dof--; nLTf--;

	/* set */	
	SetValues(theBoss, node, dof, nLTf, value);
}

void FBC_CardT::SetValues(const NodeManagerPrimitive& theBoss, int node, int dof,
	int nLTf, double value)
{
	/* set */
	fNode  = node;
	fDOF   = dof;
	fLTf   = nLTf;
	fValue = value;
	
	/* resolve the pointer to the LTf */
	fLTfPtr = theBoss.GetLTfPtr(fLTf);
}

/* return the current value */
double FBC_CardT::CurrentValue(void) const
{
	return fValue*(fLTfPtr->LoadFactor());
}

/* I/O */
void FBC_CardT::WriteHeader(ostream& out) const
{
	double* junk = NULL;
	int d_width = OutputWidth(out, junk);

	out << setw(kIntWidth) << "node" << setw(kIntWidth)    << "dof";
	out << setw(kIntWidth) << "LTf"  << setw(d_width) << "value";
	out << '\n';
}

void FBC_CardT::WriteValues(ostream& out) const
{
	double* junk = NULL;
	int d_width = OutputWidth(out, junk);

	out << setw(kIntWidth) << fNode + 1
	    << setw(kIntWidth) << fDOF + 1
<< setw(kIntWidth) << fLTf + 1
<< setw(d_width)   << fValue
<< '\n';
}
