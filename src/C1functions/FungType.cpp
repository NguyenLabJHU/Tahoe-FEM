/* $Id: FungType.cpp,v 1.5 2009-04-23 14:51:51 thao Exp $ */

#include "FungType.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

FungType::FungType(double A, double B): 
	fA(A), 
	fB(B) 
{ 
	SetName("fung-type");
}

FungType::FungType(void): 
	fA(0.0), 
	fB(0.0) 
{ 
	SetName("fung-type");
}

/* I/O */
void FungType::Print(ostream& out) const
{
	/* parameters */
	out << " fA. . . . . . . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " fB. . . . . . . . . . . . . . . . . . . . . . = " << fB << '\n';
}

void FungType::PrintName(ostream& out) const
{
	out << "    Fung Type\n";
}

/*
* Returning values
*/
double FungType::Function(double r) const
{
	return (fA/fB*(exp(fB*(r - 1.0))+ fB/r) );
}

double FungType::DFunction(double r) const
{
	return ( fA* (exp(fB*(r-1.0)) - 1.0/(r*r) ));
}

double FungType::DDFunction(double r) const
{
	return (fA*( fB*exp(fB*(r-1.0)) + 2.0/(r*r*r)) );
}

/* returning values in groups */
dArrayT& FungType::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pU++ = (fA/fB*(exp(fB*(r - 1.0))+ fB/r));
	}
	return(out);
}

dArrayT& FungType::MapDFunction(const dArrayT& in, dArrayT& out) const
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
		*pdU++ = ( fA* (exp(fB*(r-1.0)) - 1.0/(r*r) ) );
	}
	return(out);
}

dArrayT& FungType::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pddU++ = (fA*( fB*exp(fB*(r-1.0)) + 2.0/(r*r*r)) );
	}
	return(out);
}

void FungType::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	ParameterT alpha(ParameterT::Double, "alpha");
	ParameterT beta(ParameterT::Double, "beta");
	
	alpha.AddLimit(0.0,LimitT::Lower);
	beta.AddLimit(0.0,LimitT::LowerInclusive);
	list.AddParameter(alpha);
	list.AddParameter(beta);
	
	/* set the description */
	list.SetDescription("f(I) = alpha/beta*(exp(beta*(I - 1.0)) + beta/I)");	
}

void FungType::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("alpha");
	fB = list.GetParameter("beta");

}

