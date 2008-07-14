/* $Id: ElementCardT.cpp,v 1.16 2008-07-14 18:23:19 lxmota Exp $ */
/* created: paklein (05/24/1996) */
#include "ElementCardT.h"
#include <iostream.h>
#include <iomanip.h>
#include "toolboxConstants.h"
#include "dArrayT.h"

#ifndef _FRACTURE_INTERFACE_LIBRARY_
/* for the BC codes */
#include "KBC_CardT.h"
#endif

using namespace Tahoe;

/* array behavior */
namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<ElementCardT>::fByteCopy = false;
} /* namespace Tahoe */

/* initialize static data */
iArrayT ElementCardT::i_junk;
dArrayT ElementCardT::d_junk;

ElementCardT::StatusT ElementCardT::int2StatusT(int i)
{
	if (i < kOFF || i > kMarkOFF)
		ExceptionT::GeneralFail("ElementCardT::int2StatusT", "unrecognized status %d", i);
	StatusT int2status[5] = {kOFF, kON, kMarked, kMarkON, kMarkOFF};
	return int2status[i];
}

/* constructors */
ElementCardT::ElementCardT(void):
	fMatNum(-1),
	fFlag(kON),
	fNodesU(&fNodesX), // assuming isoparametric
	fData(NULL)
{

}

ElementCardT::ElementCardT(const ElementCardT& source):
	fData(NULL)
{
	/* use assignment operator */
	operator=(source);
}

/* destructor */
ElementCardT::~ElementCardT(void) { delete fData; }

/* assignment operator */
ElementCardT& ElementCardT::operator=(const ElementCardT& rhs)
{
	/* copy material number */
	fMatNum = rhs.fMatNum;
	fFlag = rhs.fFlag;

	/* shallow copies of grouped data */
	fNodesX.Alias(rhs.fNodesX);
	fEqnos.Alias(rhs.fEqnos);

	if (rhs.fNodesU == &(rhs.fNodesX))
		fNodesU = &fNodesX;    // keep isoparametric
	else
		fNodesU = rhs.fNodesU; // trust external fNodesU

	/* element storage */
	if (rhs.fData)
	{
		/* already allocated */
		if (fData)
			*fData = *(rhs.fData);
		else
		{
			fData = new ElementStorageT(*rhs.fData);
			if (!fData) throw ExceptionT::kOutOfMemory;
		}
	}

	return *this;
}

/* set material number */
void ElementCardT::SetMaterialNumber(int matnum) { fMatNum = matnum; }

/* restart operations */
void ElementCardT::ReadRestart(istream& in)
{
	int flag;
	in >> flag;
	fFlag = int2StatusT(flag);

	/* read data size */
	int i_size, d_size;
	in >> i_size >> d_size;

	/* allocate space */
	Dimension(i_size,d_size);

	/* read data */
	in >> (*fData);
}

void ElementCardT::WriteRestart(ostream& out) const
{
	out << fFlag;

	/* error to call if not allocated */
	if (!fData) throw ExceptionT::kGeneralFail;

	/* output data size */
	out << " " << IntegerData().Length();
	out << " " << DoubleData().Length();
	out << '\n';

	/* output data */
	out << (*fData);
}

/* element storage accessors/modifiers */
void ElementCardT::Dimension(int i_size, int d_size)
{
	/* nothing to do */
	if (IntegerData().Length() == i_size &&
	     DoubleData().Length() == d_size) return;

#if __option(extended_errorcheck)
	/* warning */
	if (fData != NULL) {
		cout << "\n ElementCardT::Allocate: WARNING: element data already exists\n"
		     <<   "     and will be overwritten." << endl;
	}
#endif

	/* free existing memory */
	if (fData != NULL)
	{
		delete fData;
		fData = NULL;
	}

	fData = new ElementStorageT(i_size, d_size);
	if (!fData) throw ExceptionT::kOutOfMemory;
}

void ElementCardT::Set(int i_size, int* i_data, int d_size, double* d_data)
{
	/* allocate storage card */
	if (!fData) fData = new ElementStorageT();

	/* NOTE: exisintg data gone */
	fData->Set(i_size, i_data, d_size, d_data);
}

/* make arrays alias to other data */
void ElementStorageT::Set(int i_size, int* i_data, int d_size, double* d_data)
{
	fIntegerData.Set(i_size, i_data);
	fDoubleData.Set(d_size, d_data);
}

namespace Tahoe {

/* I/O operators */
istream& operator>>(istream& in, ElementStorageT& data)
{
	in >> data.fIntegerData;
	in >> data.fDoubleData;

	return in;
}

ostream& operator<<(ostream& out, const ElementStorageT& data)
{
	out << data.fIntegerData.wrap_tight(1) << '\n'
	    << data.fDoubleData.wrap_tight(1) << '\n';

	return out;
}

}
