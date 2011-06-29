/* $Id: GaoVicky.cpp,v 1.6 2004-06-19 23:27:18 paklein Exp $ */
/* created: paklein (12/26/1998) */
#include "GaoVicky.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
GaoVicky::GaoVicky(double A, double B, double C, double D, double L):
	fA(A),
	fB(B),
	fC(C),
	fD(D),
	fL(L)
{
	SetName("Gao-Nguyen");
// should insert some consistency checks on the parameters
}

GaoVicky::GaoVicky(void):
	fA(0.0),
	fB(0.0),
	fC(0.0),
	fD(0.0),
	fL(0.0)
{
	SetName("Gao-Nguyen");
}

/* I/O */
void GaoVicky::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " << fA << '\n';
	out << "      B = " << fB << '\n';
	out << "      C = " << fC << '\n';
	out << "      D = " << fD << '\n';
	out << "      L = " << fL << '\n';
}

void GaoVicky::PrintName(ostream& out) const
{
	out << "    Gao-Nguyen\n";
}

/* returning values */
double GaoVicky::Function(double x) const
{
#pragma unused(x)
//	double dl = x - fL;
	
	cout << "\n GaoVicky::Function: only f' and f\" have been implemented\n";
	cout <<   " The function value f is not available in closed form.";
	
	throw ExceptionT::kGeneralFail;
	return 0.0;
}

double GaoVicky::DFunction(double x) const
{
	double dl = x - fL;
	double BC = fB/fC;
			
        return fA*dl/(1. + exp((-BC + dl)/fD));
}

double GaoVicky::DDFunction(double x) const
{
	double dl  =  x - fL;
	double BC = fB/fC;
	
	return fA/(1. + exp((-BC + dl)/fD)) - fA*dl*exp((-BC + dl)/fD)/
	         (fD*pow(1. + exp((-BC + dl)/fD), 2));
}

/* returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above */
dArrayT& GaoVicky::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fL;
		double BC = fB/fC;
		
		*pdU++ = fA*dl/(1. + exp((-BC + dl)/fD));
	}
	
	return out;
}

dArrayT& GaoVicky::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fL;
	        double BC = fB/fC;
		
		*pddU++ = fA/(1. + exp((-BC + dl)/fD)) - fA*dl*exp((-BC + dl)/fD)/
	         (fD*pow(1. + exp((-BC + dl)/fD), 2));
	}
	
	return out;
}

/* describe the parameters needed by the interface */
void GaoVicky::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	list.SetDescription("F(dr) = A dr/(1 + exp[(-B/C + dr)/D]");

	// should insert some consistency checks on the parameters
	list.AddParameter(fA, "A");
	list.AddParameter(fB, "B");
	list.AddParameter(fC, "C");
	list.AddParameter(fD, "D");
	list.AddParameter(fL, "L");
}

/* accept parameter list */
void GaoVicky::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("A");
	fB = list.GetParameter("B");
	fC = list.GetParameter("C");
	fD = list.GetParameter("D");
	fL = list.GetParameter("L");
}
