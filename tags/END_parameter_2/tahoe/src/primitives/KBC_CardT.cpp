/* $Id: KBC_CardT.cpp,v 1.13.6.1 2004-03-22 18:36:18 paklein Exp $ */
/* created: paklein (05/23/1996) */
#include "KBC_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "toolboxConstants.h"
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
	fvalue(0.0),
	fSchedule(NULL)
{ 

}

KBC_CardT::KBC_CardT(int node, int dof, CodeT code, const ScheduleT* schedule, double value):
	fnode(-1),
	fdof(-1),
	fcode(kFix),
	fvalue(0.0),
	fSchedule(NULL)
{
	SetValues(node, dof, code, schedule, value);
}

void KBC_CardT::SetValues(istream& in)
{
#pragma message ("delete me")
#pragma unused(in)
ExceptionT::GeneralFail();
#if 0
	/* parameters */
	int node;
	int dof;
	CodeT code;
	int nLTf;
	double value;

	/* read */
	in >> node >> dof >> code >> nLTf >> value;

	/* correct offset */
	node--; dof--; nLTf--;
	
	/* set and echo */
	SetValues(node, dof, code, nLTf, value);
#endif
}						

void KBC_CardT::SetValues(int node, int dof, CodeT code,  const ScheduleT* schedule, double value)
{
	/* set */
	fnode = node;
	fdof = dof;
	fcode = code;
	fSchedule = schedule;
	fvalue = value;

	/* fixed */
	if (fcode == kFix) {
		fSchedule = NULL;
		fvalue = 0.0;
	}
}

//void KBC_CardT::SetSchedule(const ScheduleT* schedule) { fSchedule = schedule; }

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

/* I/O */
void KBC_CardT::WriteHeader(ostream& out)
{
	int d_width = out.precision() + kDoubleExtra;

	out << setw(kIntWidth) << "node"  << setw(kIntWidth) << "dof";
//	out << setw(kIntWidth) << "code"  << setw(kIntWidth) << "LTf";
    out << setw(kIntWidth) << "code";
	out << setw(d_width)   << "value" << '\n';
}


void KBC_CardT::WriteValues(ostream& out) const
{
	int d_width = out.precision() + kDoubleExtra;

	out << setw(kIntWidth) << fnode + 1 << setw(kIntWidth) << fdof + 1;
//	out << setw(kIntWidth) << fcode     << setw(kIntWidth) << fSchedNum + 1;
	out << setw(kIntWidth) << fcode;
	out << setw(d_width)   << fvalue    << '\n';
}

KBC_CardT::CodeT KBC_CardT::int_to_CodeT (int i_code)
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

namespace Tahoe {

/* input operator for codes */
istream& operator>>(istream& in, KBC_CardT::CodeT& code)
{
	int i_code;
	in >> i_code;
	code = KBC_CardT::int_to_CodeT(i_code);
	return in;
}

}
