/* $Id: LinearT.cpp,v 1.5 2003-11-10 22:14:00 cjkimme Exp $ */
/* created: paklein (03/25/1999) */
#include "LinearT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
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
	out << "    Linear\n";
}

/* returning values in groups */
dArrayT& LinearT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

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
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

//	double* pl = in.Pointer();
//	double* pU = out.Pointer();
	
	out = fA;

    	return out;

}

dArrayT& LinearT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

//	double* pl = in.Pointer();
//	double* pU = out.Pointer();
	
	out = 0.0;

    	return out;

}
