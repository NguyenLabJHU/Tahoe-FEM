/* $Id: KBC_CardT.cpp,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
/* created: paklein (05/23/1996)                                          */

#include "KBC_CardT.h"

#include <iostream.h>
#include <iomanip.h>

#include "Constants.h"
#include "ExceptionCodes.h"

#include "LoadTime.h"

/* constructor */
KBC_CardT::KBC_CardT(void):fLTfPtr(NULL) { }
KBC_CardT::KBC_CardT(int node, int dof, CodeT code, int nLTF, double value):
	fnode(node), fdof(dof), fcode(code), fnLTf(nLTF), fvalue(value)
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
	fnode  = node;
	fdof   = dof;
	fcode  = code;
	fnLTf  = nLTf;
	fvalue = value;

	/* fixed (uncorrected) */
	if (fcode == kFix)
	{
		fnLTf  = 0;
		fvalue = 0.0;
	}
}

void KBC_CardT::SetSchedule(const LoadTime* LTfPtr) { fLTfPtr = LTfPtr; }

/* returns the value of the BC */
double KBC_CardT::Value(void) const
{
	if (fcode == kFix) /* does not have value fLTfPtr!!!! */
		return 0.0;
	else
	  {
#if __option(extended_errorcheck)
		/* double check pointers are set */
		if (!fLTfPtr) throw eGeneralFail;
#endif
		return fvalue*(fLTfPtr->LoadFactor());
	  }
}

/* I/O */
void KBC_CardT::WriteHeader(ostream& out) const
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
	out << setw(kIntWidth) << fcode     << setw(kIntWidth) << fnLTf + 1;
	out << setw(d_width)   << fvalue    << '\n';
}

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
			throw eBadInputValue;	
	}

	return in;
}
