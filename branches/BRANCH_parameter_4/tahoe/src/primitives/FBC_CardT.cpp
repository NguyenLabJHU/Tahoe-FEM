/* $Id: FBC_CardT.cpp,v 1.13.2.3 2004-07-08 07:40:24 paklein Exp $ */
/* created: paklein (06/15/1996) */
#include "FBC_CardT.h"
#include "ScheduleT.h"

#include <iostream.h>
#include <iomanip.h>

//#include "toolboxConstants.h"
//#include "NodeManagerT.h" // needed for schedule information

using namespace Tahoe;

/* copy behavior for arrays FBC_CardT's */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_CardT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_CardT>::fByteCopy = false;
} /* namespace Tahoe */

/* constructor */
FBC_CardT::FBC_CardT(void):
	fNode(-1),
	fDOF(-1),
	fValue(0.0),
	fSchedule(NULL)
{

}

#if 0
/* modifiers */
void FBC_CardT::SetValues(const NodeManagerT&, ifstreamT&)
{
#pragma message("delete me")
ExceptionT::GeneralFail();
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
#endif

void FBC_CardT::SetValues(int node, int dof, const ScheduleT* schedule, double value)
{
	/* set */
	fNode     = node;
	fDOF      = dof;
	fSchedule = schedule;
	fValue    = value;
}

/* split force value in half */
void FBC_CardT::SplitForce(void)
{
	fValue *= 0.5;
}

/* return the current value */
double FBC_CardT::CurrentValue(void) const
{
	if (!fSchedule)
		return fValue;
	else
		return fValue*(fSchedule->Value());
}

#if 0
/* I/O */
void FBC_CardT::WriteHeader(ostream& out)
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
<< setw(kIntWidth) << fSchedule
<< setw(d_width)   << fValue
<< '\n';
}
#endif
