/* $Id: FungType.cpp,v 1.9 2010-06-24 13:48:29 thao Exp $ */

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
	SetName("fung_type");
}

FungType::FungType(void): 
	fA(0.0), 
	fB(0.0) 
{ 
	SetName("fung_type");
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
	out << "   Fung Potential Squared\n";
}

/*
* Returning values
*/
double FungType::Function(double r) const
{
	/*r = I4*/
	/*Returns W(I4)*/
	double x = (r - 1.0);
	
	return (fA/(2.0*fB) * (exp(0.5*fB * x*x) - 1.0) );
}

double FungType::DFunction(double r) const
{
	/*r = I4*/
	/*Returns dW(I4)/d(I4)*/
	double x = (r - 1.0);
	
	return (0.5*fA *x * exp(0.5*fB* x*x));
}

double FungType::DDFunction(double r) const
{
	/*r = I4*/
	/*Returns dW^2(I4)/d(I4)^2*/
	double x = (r - 1.0);
	return (0.5*fA * exp(0.5*fB* x*x)*(1.0 + fB * x*x));
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
		*pU++ = Function(r);
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
		*pdU++ = DFunction(r);
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
		*pddU++ = DDFunction(r);
	}
	return(out);
}

void FungType::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	list.AddParameter(fA, "alpha");
	list.AddParameter(fB, "beta");
	
	/* set the description */
	list.SetDescription("f(I4) = alpha/beta *(exp(beta*(0.5(I - 1.0))^2) -1)");	
}

void FungType::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("alpha");
	fB = list.GetParameter("beta");

	/* check */
	if (fA < kSmall) ExceptionT::BadInputValue("ScaledSinh::TakeParameterList",
		"expecting a positive value alpha: %d", fA);
	if (fB < kSmall) ExceptionT::BadInputValue("ScaledSinh::TakeParameterList",
		"expecting a positive value for beta: %d", fB);
}

