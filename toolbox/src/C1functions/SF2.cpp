/* $Id: SF2.cpp,v 1.3 2003-11-19 22:10:39 thao Exp $ */
/* created: paklein (10/30/1997)                                          */

#include "SF2.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/*
* constructors
*/

using namespace Tahoe;

SF2::SF2(double A, double B, double l_0):
  fA(A), fB(B), fl_0(l_0) { cout << "Constructed \n";}

/*
* I/O
*/
void SF2::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " << fA << '\n';
	out << "      B = " << fB << '\n';
	out << "    l_0 = " << fl_0 << '\n';
}

void SF2::PrintName(ostream& out) const
{
	out << "    SF2\n";
}

/*
* Returning values
*/
double SF2::Function(double x) const
{
	double dl = x - fl_0;
	/*	if (dl > 0)
	  return(-(0.5*fA*fB)/exp((dl*dl)/fB));
	else
	return (0.5*fA*dl*dl-fA*fB*0.5);*/
	if (dl > 0)
	  return (-0.5*fA*fB/exp((dl*dl)/fB) + 0.5*fA*fB);
	else
	  return (0.5*fA*dl*dl);
}

double SF2::DFunction(double x) const
{
	double dl = x - fl_0;
	if (dl > 0)
	  return((fA*(dl))/exp((dl*dl)/fB));
	else
	  return(fA*dl);
}

double SF2::DDFunction(double x) const
{
	double dl = x - fl_0;
	if (dl>0)
	  return((fA*(fB - 2.0*(dl*dl)))/(fB*exp((dl*dl)/fB)));
	else
	  return(fA);
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& SF2::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = (*pl++) ;
		*pU++ = Function(x);
	}
	return(out);
}

dArrayT& SF2::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = (*pl++);
		*pdU++ = DFunction(x);
	}
	return(out);
}

dArrayT& SF2::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = (*pl++);
		*pddU++ = DDFunction(x);
	}
	return(out);
}
