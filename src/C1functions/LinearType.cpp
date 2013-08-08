/* $IdLinear: LinearType.cpp,v 1.3 2011/12/01 20:37:57 beichuan Exp $ */

#include "LinearType.h"
#include <iostream>
#include <cmath>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

LinearType::LinearType(double A): 
	fA(A)
{ 
	SetName("linear_type");
}

LinearType::LinearType(void): 
	fA(0.0)
{ 
	SetName("linear_type");
}

/* I/O */
void LinearType::Print(ostream& out) const
{
	/* parameters */
	out << " fA. . . . . . . . . . . . . . . . . . . . . = " << fA << '\n';
}

void LinearType::PrintName(ostream& out) const
{
	out << "    linear type \n";
}

/*
* Returning values
*/
double LinearType::Function(double r) const
{
	return (fA/2*(r-1.0)*(r-1.0));
}

double LinearType::DFunction(double r) const
{
	return (fA*(r-1.0));
}

double LinearType::DDFunction(double r) const
{
	return (fA);
}

/* returning values in groups */
dArrayT& LinearType::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pU++ = (fA/2*(r-1.0)*(r-1.0));
	}
	return(out);
}

dArrayT& LinearType::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
//		cout <<"\nr: "<<r;
//		cout<<"\ndU: "<<(r-1.0);
		*pdU++ = (fA*(r-1.0));
	}
	return(out);
}

dArrayT& LinearType::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pddU++ = (fA);
	}
	return(out);
}

void LinearType::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	ParameterT theta(ParameterT::Double, "theta");
	
	
	theta.AddLimit(0.0,LimitT::Lower);
	list.AddParameter(theta);
	
	/* set the description */
	list.SetDescription("f(I) = theta/2*(I4-1)^2");	
}

void LinearType::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("theta");

}

