/* $Id: ElementCardT.cpp,v 1.4 2001-09-05 21:56:41 paklein Exp $ */
/* created: paklein (05/24/1996)                                          */

#include "ElementCardT.h"
#include <iostream.h>
#include <iomanip.h>
#include "Constants.h"
#include "dArrayT.h"

/* for the BC codes */
#include "KBC_CardT.h"

/* array behavior */
const bool ArrayT<ElementCardT>::fByteCopy = false;

/* initialize static data */
iArrayT ElementCardT::i_junk;
dArrayT ElementCardT::d_junk;

/* constructors */
ElementCardT::ElementCardT(void):
	fMatNum(-1),
	fFlag(1),
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
			if (!fData) throw(eOutOfMemory);
		}
	}

	return *this;
}

/* set material number */
void ElementCardT::SetMaterialNumber(int matnum) { fMatNum = matnum; }

/* restart operations */
void ElementCardT::ReadRestart(istream& in)
{
	in >> fFlag;

	/* read data size */
	int i_size, d_size;
	in >> i_size >> d_size;
	
	/* allocate space */
	Allocate(i_size,d_size);

	/* read data */
	in >> (*fData);
}

void ElementCardT::WriteRestart(ostream& out) const
{
	out << fFlag;

	/* error to call if not allocated */
	if (!fData) throw eGeneralFail;

	/* output data size */
	out << " " << IntegerData().Length();
	out << " " << DoubleData().Length();
	out << '\n';

	/* output data */
	out << (*fData);
}

/* element storage accessors/modifiers */
void ElementCardT::Allocate(int i_size, int d_size)
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
	if (!fData) throw eOutOfMemory;
}

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
