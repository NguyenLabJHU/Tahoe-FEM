/* $Id: SmithFerrante.cpp,v 1.6 2004-03-24 01:56:34 paklein Exp $ */
/* created: paklein (10/30/1997) */
#include "SmithFerrante.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
SmithFerrante::SmithFerrante(double A, double B, double l_0):
	fA(A), 
	fB(B), 
	fl_0(l_0) 
{ 
	SetName("Smith-Ferrante");
}

SmithFerrante::SmithFerrante(void):
	fA(0), 
	fB(0), 
	fl_0(0) 
{ 
	SetName("Smith-Ferrante");
}

/* I/O */
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

/* returning values */
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

/* returning values in groups */
dArrayT& SmithFerrante::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
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
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
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
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fl_0;
		*pddU++ = (fA*(fB - (dl)))/(fB*exp((dl)/fB));
	}
	return(out);
}

/* describe the parameters needed by the interface */
void SmithFerrante::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	list.AddParameter(fA, "a");
	list.AddParameter(fB, "b");
	list.AddParameter(fl_0, "x_0");

	/* set the description */
	list.SetDescription("f(x) = a*(x - x_0)*exp(-(x - x_0)/B)");
}

/* accept parameter list */
void SmithFerrante::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);
	
	fA = list.GetParameter("a");
	fB = list.GetParameter("b");
	fl_0 = list.GetParameter("x_0");

}
