/* $Id: VWType.cpp,v 1.2 2010-06-24 13:48:30 thao Exp $ */

#include "VWType.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

VWType::VWType(double A, double B): 
	fA(A), 
	fB(B) 
{ 
	SetName("vw-type");
}

VWType::VWType(void): 
	fA(0.0), 
	fB(0.0) 
{ 
	SetName("vw-type");
}

/* I/O */
void VWType::Print(ostream& out) const
{
	/* parameters */
	out << " fA. . . . . . . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " fB. . . . . . . . . . . . . . . . . . . . . . = " << fB << '\n';
}

void VWType::PrintName(ostream& out) const
{
	out << "    Veronda Westmann type Type\n";
}

/*
* Returning values
*/
double VWType::Function(double r) const
{
	return (fA/fB*(exp(fB*(r - 1.0))- fB*r) );
}

double VWType::DFunction(double r) const
{
	return ( fA* (exp(fB*(r-1.0)) - 1.0 ));
}

double VWType::DDFunction(double r) const
{
	return (fA*( fB*exp(fB*(r-1.0)) ) );
}

/* returning values in groups */
dArrayT& VWType::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pU++ = (fA/fB*(exp(fB*(r - 1.0))- fB*r));
	}
	return(out);
}

dArrayT& VWType::MapDFunction(const dArrayT& in, dArrayT& out) const
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
		*pdU++ = ( fA* (exp(fB*(r-1.0)) - 1.0 ) );
	}
	return(out);
}

dArrayT& VWType::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pddU++ = (fA*( fB*exp(fB*(r-1.0)) ) );
	}
	return(out);
}

void VWType::DefineParameters(ParameterListT& list) const
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
	list.SetDescription("f(I) = alpha/beta*(exp(beta*(I - 1.0)) + beta*I)");	
}

void VWType::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("alpha");
	fB = list.GetParameter("beta");

}

