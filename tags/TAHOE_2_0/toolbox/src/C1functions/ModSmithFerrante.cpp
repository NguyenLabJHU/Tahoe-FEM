/* $Id: ModSmithFerrante.cpp,v 1.6 2004-04-27 07:22:19 paklein Exp $ */

/* Smith Ferrante modified to have a linear branch */

#include "ModSmithFerrante.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
ModSmithFerrante::ModSmithFerrante(void):
	fA(0.0),
	fB(0.0) 
{ 
	SetName("modified_Smith-Ferrante");
}

ModSmithFerrante::ModSmithFerrante(double A, double B):
	fA(A), 
	fB(B) 
{ 
	SetName("modified_Smith-Ferrante");
}

/* I/O */
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

/* returning values */
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

	const double* pl = in.Pointer();
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

	const double* pl  = in.Pointer();
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

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	double x;
	
	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++;
		*pddU++ = (x > 0) ? ((fA*(x))/exp((x)/fB)) : fA*x;
	}
	return(out);
}

/* describe the parameters needed by the interface */
void ModSmithFerrante::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	list.AddParameter(fA, "a");
	list.AddParameter(fB, "b");
}

/* accept parameter list */
void ModSmithFerrante::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);
	fA = list.GetParameter("a");
	fB = list.GetParameter("b");
}
