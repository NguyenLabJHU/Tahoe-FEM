/* $Id: ModSmithFerrante.cpp,v 1.4 2003-11-11 01:49:32 rjones Exp $ */

/* Smith Ferrante modified to have a linear branch */

#include "ModSmithFerrante.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/*
* constructors
*/

using namespace Tahoe;

ModSmithFerrante::ModSmithFerrante(double A, double B):
	fA(A), fB(B) { }

/*
* I/O
*/
void ModSmithFerrante::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " << fA << '\n';
	out << "      B = " << fB << '\n';
}

void ModSmithFerrante::PrintName(ostream& out) const
{
	out << "    Smith-Ferrante\n";
}

/*
* Returning values
*/
double ModSmithFerrante::Function(double x) const
{
	return (x > 0) ? (-((fA*fB*(fB + (x)))/exp((x)/fB))) : 0.5*fA*x*x;
}

double ModSmithFerrante::DFunction(double x) const
{
	return (x > 0) ? ((fA*(x))/exp((x)/fB)) : fA*x;
}

double ModSmithFerrante::DDFunction(double x) const
{
	return (x > 0) ? ((fA*(fB - (x)))/(fB*exp((x)/fB))) : fA;
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& ModSmithFerrante::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x;
	
	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pU++ = (x > 0) ? (-((fA*fB*(fB + (x)))/exp((x)/fB))) : 0.5*fA*x*x;
	}
	return(out);
}

dArrayT& ModSmithFerrante::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	double x;
	
	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pdU++ = (x > 0) ? ((fA*(x))/exp((x)/fB)) : fA*x;
	}
	return(out);
}

dArrayT& ModSmithFerrante::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	double x;
	
	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pddU++ = (x > 0) ? ((fA*(x))/exp((x)/fB)) : fA*x;
	}
	return(out);
}
