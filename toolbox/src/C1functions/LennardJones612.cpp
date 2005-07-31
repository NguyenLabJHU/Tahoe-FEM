/* $Id: LennardJones612.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (10/30/1997)                                          */

#include "LennardJones612.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"

/* constants */
const double twoe1by6 = pow(2.0,1.0/6.0);

/*
* constructors
*/
LennardJones612::LennardJones612(double A): fA(A) { }

/*
* I/O
*/
void LennardJones612::Print(ostream& out) const
{
	/* parameters */
	out << " Scaling constant. . . . . . . . . . . . . . . . = " << fA << '\n';
}

void LennardJones612::PrintName(ostream& out) const
{
	out << "    Lennard-Jones 6-12\n";
}

/*
* Returning values
*/
double LennardJones612::Function(double x) const
{
	return fA*(0.5*pow(x,-12.0) - pow(x,-6.0));
}

double LennardJones612::DFunction(double x) const
{
	return 6.0*fA*(-pow(x,-13.0) + pow(x,-7.0));
}

double LennardJones612::DDFunction(double x) const
{
	return fA*(78.0*pow(x,-14.0) - 42.0*pow(x,-8.0));
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& LennardJones612::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;					
		*pU++ = fA*(0.5*pow(r,-12.0) - pow(r,-6.0));
	}
	return out;
}

dArrayT& LennardJones612::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;					
		*pdU++ = 6.0*fA*(-pow(r,-13.0) + pow(r,-7.0));
	}
	return out;
}

dArrayT& LennardJones612::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;				
		*pddU++ = fA*(78.0*pow(r,-14.0) - 42.0*pow(r,-8.0));
	}
	return out;
}
