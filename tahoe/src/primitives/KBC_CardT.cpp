/* $Id: KBC_CardT.cpp,v 1.14.30.1 2005-07-25 02:37:24 paklein Exp $ */
/* created: paklein (05/23/1996) */
#include "KBC_CardT.h"
#include "ExceptionT.h"
#include "ScheduleT.h"

using namespace Tahoe;

/* copy behavior for arrays KBC_CardT's */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<KBC_CardT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<KBC_CardT>::fByteCopy = false;
} /* namespace Tahoe */

/* constructor */
KBC_CardT::KBC_CardT(void):
	fnode(-1),
	fdof(-1),
	fcode(kFix),
	fmode(kUndefined),
	fvalue(0.0),
	fSchedule(NULL)
{ 

}

#pragma message("delete me")
#if 0
KBC_CardT::KBC_CardT(int node, int dof, CodeT code, const ScheduleT* schedule, double value):
	fnode(-1),
	fdof(-1),
	fcode(kFix),
	fmode(kNode),	
	fvalue(0.0),
	fSchedule(NULL)
{
	SetValues(node, dof, code, schedule, value);
}
#endif

void KBC_CardT::SetValues(int node, int dof, CodeT code,  const ScheduleT* schedule, double value)
{
	/* set */
	fnode = node; 
	fID.Clear();
	fdof = dof;
	fcode = code;
	fSchedule = schedule;
	fvalue = value;
	fmode= kNode;

	/* fixed */
	if (fcode == kFix) {
		fSchedule = NULL;
		fvalue = 0.0;
	}
}

void KBC_CardT::SetValues(const StringT& ID, int dof, CodeT code,  const ScheduleT* schedule, double value)
{
	/* set */
	fnode = -1;
	fID = ID;
	fdof = dof;
	fcode = code;
	fSchedule = schedule;
	fvalue = value;
	fmode= kSet;

	/* fixed */
	if (fcode == kFix) {
		fSchedule = NULL;
		fvalue = 0.0;
	}
}

/* returns the value of the BC */
double KBC_CardT::Value(void) const
{
	if (fcode == kFix) /* does not have value fLTfPtr!!!! */
		return 0.0;
	else
	{
		if (!fSchedule)
			return fvalue;
		else
			return fvalue*(fSchedule->Value());
	  }
}

KBC_CardT::CodeT KBC_CardT::int2CodeT(int i_code)
{
	/* resolve code */
	switch (i_code)
	{
		case KBC_CardT::kFix:
			return KBC_CardT::kFix;
			break;
		case KBC_CardT::kDsp:
			return KBC_CardT::kDsp;
			break;
		case KBC_CardT::kVel:
			return KBC_CardT::kVel;
			break;
		case KBC_CardT::kAcc:
			return KBC_CardT::kAcc;
			break;
		default:
			ExceptionT::BadInputValue("KBC_CardT::int_to_CodeT", "unknown code: %d", i_code);
	}

	/* dummy */
	return KBC_CardT::kAcc;
}
