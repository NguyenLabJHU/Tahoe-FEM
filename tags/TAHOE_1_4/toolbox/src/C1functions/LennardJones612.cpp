/* $Id: LennardJones612.cpp,v 1.7 2004-05-10 17:01:28 paklein Exp $ */
/* created: paklein (10/30/1997) */
#include "LennardJones612.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constants */
const double twoe1by6 = pow(2.0,1.0/6.0);

/* constructors */
LennardJones612::LennardJones612(void): 
	fA(0.0),
	fB(1.0)
{ 
	SetName("Lennard-Jones_6-12");
}

LennardJones612::LennardJones612(double A, double B): 
	fA(A),
	fB(B)
{
	SetName("Lennard-Jones_6-12");
}

/* I/O */
void LennardJones612::Print(ostream& out) const
{
	/* parameters */
	out << " Scaling constant. . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " Length scaling. . . . . . . . . . . . . . . . . = " << fB << '\n';
}

void LennardJones612::PrintName(ostream& out) const
{
	out << "    Lennard-Jones 6-12\n";
}

/* returning values */
double LennardJones612::Function(double x) const
{
	double s = x/fB;
	return fA*(0.5*pow(s,-12.0) - pow(s,-6.0));
}

double LennardJones612::DFunction(double x) const
{
	double s = x/fB;
	return 6.0*fA*(-pow(s,-13.0) + pow(s,-7.0));
}

double LennardJones612::DDFunction(double x) const
{
	double s = x/fB;
	return fA*(78.0*pow(s,-14.0) - 42.0*pow(s,-8.0));
}

/* returning values in groups */
dArrayT& LennardJones612::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("LennardJones612::MapFunction");

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++/fB;
		*pU++ = fA*(0.5*pow(r,-12.0) - pow(r,-6.0));
	}
	return out;
}

dArrayT& LennardJones612::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("LennardJones612::MapDFunction");

	const double* pl = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++/fB;
		*pdU++ = 6.0*fA*(-pow(r,-13.0) + pow(r,-7.0));
	}
	return out;
}

dArrayT& LennardJones612::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) ExceptionT::GeneralFail("LennardJones612::MapDDFunction");

	const double* pl = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++/fB;
		*pddU++ = fA*(78.0*pow(r,-14.0) - 42.0*pow(r,-8.0));
	}
	return out;
}

void LennardJones612::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fA, "a");
	list.AddParameter(fB, "b");
	
	/* set the description */
	list.SetDescription("f(x) = a*((b/x)^12 - (b/x)^6)");	
}

void LennardJones612::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("a");
	fB = list.GetParameter("b");	
}
