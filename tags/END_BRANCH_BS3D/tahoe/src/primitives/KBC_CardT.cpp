/* $Id: KBC_CardT.cpp,v 1.13 2003-11-04 01:32:06 paklein Exp $ */
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
KBC_CardT::KBC_CardT(void):fSchedule(NULL) { }
KBC_CardT::KBC_CardT(int node, int dof, CodeT code, int nLTF, double value):
	fnode(node), fdof(dof), fcode(code), fSchedNum(nLTF), fvalue(value)
{

}

void KBC_CardT::SetValues(istream& in)
{
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
}						

void KBC_CardT::SetValues(int node, int dof, CodeT code, int nLTf, double value)
{
	/* set */
	fnode = node;
	fdof = dof;
	fcode = code;
	fSchedNum = nLTf;
	fvalue = value;

	/* fixed (uncorrected) */
	if (fcode == kFix)
	{
		fSchedNum  = 0;
		fvalue = 0.0;
	}
}

void KBC_CardT::SetSchedule(const ScheduleT* schedule) { fSchedule = schedule; }

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
	out << setw(kIntWidth) << "code"  << setw(kIntWidth) << "LTf";
	out << setw(d_width)   << "value" << '\n';
}


void KBC_CardT::WriteValues(ostream& out) const
{
	int d_width = out.precision() + kDoubleExtra;

	out << setw(kIntWidth) << fnode + 1 << setw(kIntWidth) << fdof + 1;
	out << setw(kIntWidth) << fcode     << setw(kIntWidth) << fSchedNum + 1;
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
			cout << "\n KBC_CardT::int_to_CodeT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
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

	/* resolve code */
	switch (i_code)
	{
		case KBC_CardT::kFix:
			code = KBC_CardT::kFix;
			break;
		case KBC_CardT::kDsp:
			code = KBC_CardT::kDsp;
			break;
		case KBC_CardT::kVel:
			code = KBC_CardT::kVel;
			break;
		case KBC_CardT::kAcc:
			code = KBC_CardT::kAcc;
			break;
		default:
			cout << "\n operator>>KBC_CardT::CodeT: unknown code: "
			<< i_code<< endl;
			throw ExceptionT::kBadInputValue;	
	}

	return in;
}

}
