/* $Id: LinearT.cpp,v 1.2 2002-07-02 19:56:31 cjkimme Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "LinearT.h"
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"

/* constructors */

using namespace Tahoe;

LinearT::LinearT(double A, double B): 
  fA(A),
  fB(B)
{ }

/* I/O */
void LinearT::Print(ostream& out) const
{
	/* parameters */
	out <<"A: "<< fA << '\n';
	out <<"B: "<< fB << '\n';
}

void LinearT::PrintName(ostream& out) const
{
        out << "Function .............  linear "<<'\n';
}

/* returning values in groups */
dArrayT& LinearT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x;
	
	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++; 
		*pU++ = Function(x);
	}

	return out;
}

dArrayT& LinearT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	out = fA;

    	return out;

}

dArrayT& LinearT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	out = 0.0;

    	return out;

}
