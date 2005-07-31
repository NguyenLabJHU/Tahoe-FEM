/* $Id: SmithFerrante.cpp,v 1.1.1.1 2001-01-25 20:56:27 paklein Exp $ */
/* created: paklein (10/30/1997)                                          */

#include "SmithFerrante.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"

/*
* constructors
*/
SmithFerrante::SmithFerrante(double A, double B, double l_0):
	fA(A), fB(B), fl_0(l_0) { }

/*
* I/O
*/
void SmithFerrante::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " << fA << '\n';
	out << "      B = " << fB << '\n';
	out << "    l_0 = " << fl_0 << '\n';
}

void SmithFerrante::PrintName(ostream& out) const
{
	out << "    Smith-Ferrante\n";
}

/*
* Returning values
*/
double SmithFerrante::Function(double x) const
{
	double dl = x - fl_0;
	return(-((fA*fB*(fB + (dl)))/exp((dl)/fB)));
}

double SmithFerrante::DFunction(double x) const
{
	double dl = x - fl_0;
	return((fA*(dl))/exp((dl)/fB));
}

double SmithFerrante::DDFunction(double x) const
{
	double dl = x - fl_0;
	return((fA*(fB - (dl)))/(fB*exp((dl)/fB)));
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& SmithFerrante::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fl_0;
		*pU++ = -((fA*fB*(fB + (dl)))/exp((dl)/fB));
	}
	return(out);
}

dArrayT& SmithFerrante::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fl_0;
		*pdU++ = (fA*(dl))/exp((dl)/fB);
	}
	return(out);
}

dArrayT& SmithFerrante::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fl_0;
		*pddU++ = (fA*(fB - (dl)))/(fB*exp((dl)/fB));
	}
	return(out);
}
