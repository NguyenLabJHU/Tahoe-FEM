/* $Id: GaoJi2.cpp,v 1.5 2004-06-19 23:27:18 paklein Exp $ */
/* created: Baohua Ji (02/25/2002) */
#include "GaoJi2.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
GaoJi2::GaoJi2(double A, double B, double C, double L_0):
	fA(A),
	fB(B),
	fN(C),
	fL_0(L_0)
{
	SetName("Gao-Ji_2");
	// should insert some consistency checks on the parameters
}

GaoJi2::GaoJi2(void):
	fA(0.0),
	fB(0.0),
	fN(0.0),
	fL_0(0.0)
{
	SetName("Gao-Ji_2");
}

/* I/O */
void GaoJi2::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " <<   fA << '\n';
	out << "      B = " <<   fB << '\n';
	out << "      C = " <<   fN << '\n';
	out << "    L_0 = " << fL_0 << '\n';
}

void GaoJi2::PrintName(ostream& out) const
{
	out << "    Gao-Ji2\n";
}

/* returning values */
double GaoJi2::Function(double x) const
{
	double dl = x - fL_0;

          double B1 = fB*(pow(2.*fN,fN)*pow(-1. + 2.*fN,0.5 - fN))/2.718281828;
        if(fN > 1.)


          return( -(fA*pow(B1,2.)*pow(1. + pow(dl,2.)/pow(B1,2.),1. - fN))/(2.*(-1. + fN)) );

        else
          return 0.0;

}

double GaoJi2::DFunction(double x) const
{
	double dl  =  x - fL_0;
	
          double B1 = fB*(pow(2.*fN,fN)*pow(-1. + 2.*fN,0.5 - fN))/2.718281828;
        if(fN > 1.)


          return( (fA*dl)/pow(1. + pow(dl,2.)/pow(B1,2.),fN) );

        else
          return 0.0;

}

double GaoJi2::DDFunction(double x) const
{
	double dl  =  x - fL_0;

          double B1 = fB*(pow(2.*fN,fN)*pow(-1. + 2.*fN,0.5 - fN))/2.718281828;
	if(fN > 1.)


          return ( fA/pow(1. + pow(dl,2.)/pow(B1,2.),fN) - 
		   (2.*fA*pow(dl,2.)*pow(1. + pow(dl,2.)/pow(B1,2.),-1. - fN)*fN)/pow(B1,2.) );

        else
          return 0.0;

}

/* returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above */

dArrayT& GaoJi2::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();

        double B1 = fB*(pow(2.*fN,fN)*pow(-1. + 2.*fN,0.5 - fN))/2.718281828;
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl = (*pl++) - fL_0;

                
		*pU++ = (fN > 1.) ? -(fA*pow(B1,2.)*pow(1. + pow(dl,2.)/pow(B1,2.),1. - fN))/
		                         (2.*(-1. + fN)) : 0.0;
	}
	return(out);
}

dArrayT& GaoJi2::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();

        double B1 = fB*(pow(2.*fN,fN)*pow(-1. + 2.*fN,0.5 - fN))/2.718281828;
	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl  = (*pl++) - fL_0;

		
		*pdU++ = (fN > 1.) ? (fA*dl)/pow(1. + pow(dl,2.)/pow(B1,2.),fN) : 0.0;
	}
	
	return out;
}

dArrayT& GaoJi2::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();

        double B1 = fB*(pow(2.*fN,fN)*pow(-1. + 2.*fN,0.5 - fN))/2.718281828;

	
	for (int i = 0; i < in.Length(); i++)
	{
		double dl  = (*pl++) - fL_0;
		
		*pddU++ = (fN > 1.) ? fA/pow(1. + pow(dl,2.)/pow(B1,2.),fN) - 
		   (2.*fA*pow(dl,2.)*pow(1. + pow(dl,2.)/pow(B1,2.),-1. - fN)*fN)/pow(B1,2.) : 0.0;
	}
	
	return out;
}

/* describe the parameters needed by the interface */
void GaoJi2::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	// should insert some consistency checks on the parameters
	list.AddParameter(fA, "A");
	list.AddParameter(fB, "B");
	list.AddParameter(fN, "N");
	list.AddParameter(fL_0, "L_0");
}

/* accept parameter list */
void GaoJi2::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("A");
	fB = list.GetParameter("B");
	fN = list.GetParameter("N");
	fL_0 = list.GetParameter("L_0");
}
