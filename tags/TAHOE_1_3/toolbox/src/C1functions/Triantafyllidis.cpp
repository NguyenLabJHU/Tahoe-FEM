/* $Id: Triantafyllidis.cpp,v 1.4 2003-11-21 22:41:27 paklein Exp $ */

#include "Triantafyllidis.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructor */

using namespace Tahoe;

Triantafyllidis::Triantafyllidis(double A): fA(A) { }

/* I/O */
void Triantafyllidis::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " << fA << '\n';
}

void Triantafyllidis::PrintName(ostream& out) const
{
	out << "    Triantafyllidis\n";
}

/* returning values */
double Triantafyllidis::Function(double x) const
{
	if (x < 0.0) throw ExceptionT::kGeneralFail;
	return fA*(log(x) - (x - 1.0)/x);
}

double Triantafyllidis::DFunction(double x) const
{
	return fA*(x - 1.0)/x/x;
}

double Triantafyllidis::DDFunction(double x) const
{
	return fA*(2.0 - x)/x/x/x;
}

/* returning values in groups */
dArrayT& Triantafyllidis::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pl++;
		if (x < 0.0) throw ExceptionT::kGeneralFail;
		*pU++ = fA*(log(x) - (x - 1.0)/x);
	}
	return out;
}

dArrayT& Triantafyllidis::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pl++;
		*pdU++ = fA*(x - 1.0)/x/x;
	}
	return out;
}

dArrayT& Triantafyllidis::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pl++;
		*pddU++ = fA*(2.0 - x)/x/x/x;
	}
	return out;
}
