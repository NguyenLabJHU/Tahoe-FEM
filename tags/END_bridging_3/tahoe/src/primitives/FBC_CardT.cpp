/* $Id: FBC_CardT.cpp,v 1.10.20.1 2003-10-27 22:22:52 paklein Exp $ */
/* created: paklein (06/15/1996) */
#include "FBC_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
#include "NodeManagerT.h" // needed for schedule information
#include "fstreamT.h"
#include "ScheduleT.h"

using namespace Tahoe;

/* copy behavior for arrays FBC_CardT's */
namespace Tahoe {
template<> const bool ArrayT<FBC_CardT*>::fByteCopy = true;
template<> const bool ArrayT<FBC_CardT>::fByteCopy = false;
} /* namespace Tahoe */

/* constructor */
FBC_CardT::FBC_CardT(void):
	fNode(-1),
	fDOF(-1),
	fSchedNum(-1),
	fValue(0.0),
	fSchedule(NULL)
{

}

/* modifiers */
void FBC_CardT::SetValues(const NodeManagerT& theBoss, ifstreamT& in)
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

void FBC_CardT::SetValues(const NodeManagerT& theBoss, int node, int dof,
	int schedule, double value)
{
	/* set */
	fNode     = node;
	fDOF      = dof;
	fSchedNum = schedule;
	fValue    = value;
	
	/* resolve the pointer to the schedule */
	fSchedule = theBoss.Schedule(fSchedNum);
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
<< setw(kIntWidth) << fSchedNum + 1
<< setw(d_width)   << fValue
<< '\n';
}
