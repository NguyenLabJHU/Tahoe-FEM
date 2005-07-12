/* $Id: FBC_CardT.cpp,v 1.14.20.1 2005-06-08 17:21:42 paklein Exp $ */
/* created: paklein (06/15/1996) */
#include "FBC_CardT.h"
#include "ScheduleT.h"

using namespace Tahoe;

/* copy behavior for arrays FBC_CardT's */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_CardT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<FBC_CardT>::fByteCopy = false;
} /* namespace Tahoe */

/* constructor */
FBC_CardT::FBC_CardT(void):
	fNode(-1),
	fMode(kUndefined),
	fDOF(-1),
	fValue(0.0),
	fSchedule(NULL)
{

}

void FBC_CardT::SetValues(int node, int dof, const ScheduleT* schedule, double value)
{
	/* set */
	fNode     = node;
	fDOF      = dof;
	fSchedule = schedule;
	fValue    = value;
	fMode     = kNode;
	fID.Clear();
}

void FBC_CardT::SetValues(const StringT& ID, int dof, const ScheduleT* schedule, double value)
{
	/* set */
	fNode     = -1;
	fDOF      = dof;
	fSchedule = schedule;
	fValue    = value;
	fMode     = kSet;
	fID       = ID;
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
