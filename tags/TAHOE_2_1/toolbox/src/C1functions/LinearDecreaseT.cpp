/* $Id: LinearDecreaseT.cpp,v 1.5 2004-03-16 05:37:12 paklein Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "LinearDecreaseT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructors */

using namespace Tahoe;

LinearDecreaseT::LinearDecreaseT(double A, double L): 
  fA(A),
  fL(L)
{ }

/* I/O */
void LinearDecreaseT::Print(ostream& out) const
{
	/* parameters */
	out <<"A: "<< fA << " ...... L: "<<fL<< '\n';
}

void LinearDecreaseT::PrintName(ostream& out) const
{
        out << "Function .............  1-A(x/L)"<<'\n';
}

/* returning values in groups */
dArrayT& LinearDecreaseT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		*pU++ = (1-fA*(*pl)/fL);
		pl++;
	}

	return out;
}

dArrayT& LinearDecreaseT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	out = -fA/fL;

	return out;
}

dArrayT& LinearDecreaseT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	out = 0.0;
	return out;
}
