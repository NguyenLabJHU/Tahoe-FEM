/* $Id: GaoKlein.cpp,v 1.5 2003-11-21 22:41:27 paklein Exp $ */
/* created: paklein (12/26/1998)                                          */

#include "GaoKlein.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructor */

using namespace Tahoe;

GaoKlein::GaoKlein(double A, double B, double C, double L):
	fA(A),
	fB(B),
	fC(C),
	fL(L)
{
// should insert some consistency checks on the parameters
}

/* I/O */
void GaoKlein::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " << fA << '\n';
	out << "      B = " << fB << '\n';
	out << "      C = " << fC << '\n';
	out << "      L = " << fL << '\n';
}

void GaoKlein::PrintName(ostream& out) const
{
	out << "    Gao-Klein\n";
}

/* returning values */
double GaoKlein::Function(double x) const
{
#pragma unused(x)
//	double dr = x - fL;
	
	cout << "\n GaoKlein::Function: only f' and f\" have been implemented\n";
	cout <<   " The function value f is not available in closed form, but is\n";
	cout <<   " given by:\n";
cout << "                                          2";
cout << "                                       A l    l (A B - 2 A L)";
cout << "     2          2      2          2    ---- + ---------------";
cout << "-(A B  C + A B C  + A B  L - A B L )    2            2";
cout << "------------------------------------ + ---------------------- + ";
cout << "              B/(C - l + L)                 B/(C - l + L)";
cout << "         2 B E                             E";
cout << " ";
cout << "        2                               B";
cout << "  (-(A B ) - 2 A B C) ExpIntegralEi[----------]";
cout << "                                    -C + l - L";
cout << "  ---------------------------------------------";
cout << "                        2	";
cout << endl;
	
	throw ExceptionT::kGeneralFail;
	return 0.0;
}

double GaoKlein::DFunction(double x) const
{
	double dr  =  x - fL;
	double cut = dr - fC;
	
	if (cut < 0.0)
		return fA*dr*exp(fB/cut);
	else
		return 0.0;
}

double GaoKlein::DDFunction(double x) const
{
	double dr  =  x - fL;
	double cut = dr - fC;
	
	if (cut < 0.0)
		return fA*exp(fB/cut)*(1.0 - dr*fB/(cut*cut));
	else
		return 0.0;
}

/* returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above */
dArrayT& GaoKlein::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dr  = (*pl++) - fL;
		double cut = dr - fC;
		
		*pdU++ = (cut < 0.0) ? fA*dr*exp(fB/cut) : 0.0;
	}
	
	return out;
}

dArrayT& GaoKlein::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dr  = (*pl++) - fL;
		double cut = dr - fC;
		
		*pddU++ = (cut < 0.0) ? fA*exp(fB/cut)*(1.0 - dr*fB/(cut*cut)) : 0.0;
	}
	
	return out;
}
