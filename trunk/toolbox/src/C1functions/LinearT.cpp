/* $Id: LinearT.cpp,v 1.9 2004-06-09 06:23:51 paklein Exp $ */
/* created: paklein (03/25/1999) */
#include "LinearT.h"
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
LinearT::LinearT(double A, double B): 
  fA(A),
  fB(B)
{
	SetName("linear_function");
}

LinearT::LinearT(void): 
  fA(0.0),
  fB(0.0)
{ 
	SetName("linear_function");
}

/* I/O */
void LinearT::Print(ostream& out) const
{
	/* parameters */
	out <<"A: "<< fA << '\n';
	out <<"B: "<< fB << '\n';
}

void LinearT::PrintName(ostream& out) const
{
	out << "    Linear\n";
}

/* returning values in groups */
dArrayT& LinearT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x;
	
	for (int i = 0; i < in.Length(); i++)
	{
		x = *pl++; 
		*pU++ = Function(x);
	}

	return out;
}

dArrayT& LinearT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;
	out = fA;
	return out;
}

dArrayT& LinearT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;
	out = 0.0;
	return out;
}

void LinearT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fA, "a");
	list.AddParameter(fB, "b");
	
	/* set the description */
	list.SetDescription("f(x) = a*x + b");	
}

void LinearT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("a");
	fB = list.GetParameter("b");
}
