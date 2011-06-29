/* $Id: SF2.cpp,v 1.5 2004-06-19 23:27:18 paklein Exp $ */
/* created: paklein (10/30/1997) */
#include "SF2.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructors */
SF2::SF2(double A, double B, double l_0):
	fA(A), 
	fB(B), 
	fl_0(l_0)
{
	SetName("Smith-Ferrante_2");
}

SF2::SF2(void):
	fA(0.0), 
	fB(0.0), 
	fl_0(0.0)
{
	SetName("Smith-Ferrante_2");
}

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

	const double* pl = in.Pointer();
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

	const double* pl  = in.Pointer();
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

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = (*pl++);
		*pddU++ = DDFunction(x);
	}
	return(out);
}

/* describe the parameters needed by the interface */
void SF2::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);
	
	list.SetDescription("F(dr) = A dr exp(-dr^2/B)");

	// should insert some consistency checks on the parameters
	list.AddParameter(fA, "A");
	list.AddParameter(fB, "B");
	list.AddParameter(fl_0, "l_0");
}

/* accept parameter list */
void SF2::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("A");
	fB = list.GetParameter("B");
	fl_0 = list.GetParameter("l_0");
}
