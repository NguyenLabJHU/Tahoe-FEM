/* $Id: ParabolaPotT.cpp,v 1.1 2003-05-21 16:04:08 thao Exp $ */
/* created: paklein (03/25/1999)                                          */

#include "ParabolaPotT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructors */

using namespace Tahoe;

ParabolaPotT::ParabolaPotT(double k, double B, double l0): fk(k),fl0(l0) { }

/* I/O */
void ParabolaPotT::Print(ostream& out) const
{
	/* parameters */
	out << " U'' . . . . . . . . . . . . . . . . . . . . . . = " << fk << '\n';
}

void ParabolaPotT::PrintName(ostream& out) const
{
	out << "    Quadratic function\n";
}

/* returning values in groups */
dArrayT& ParabolaPotT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pl-fl0;
		*pU++ = 0.5*fk*x*x-0.5*fk*fB;
		pl++;
	}

	return out;
}

dArrayT& ParabolaPotT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pl - fl0;
		*pdU++ = fk*x;
		pl++;
	}
	return out;
}

dArrayT& ParabolaPotT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	out = fk;
	return out;
}
