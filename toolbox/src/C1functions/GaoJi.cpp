/* $Id: GaoJi.cpp,v 1.5 2004-06-19 23:27:18 paklein Exp $ */
/* created: Baohua Ji (02/25/2002) */
#include "GaoJi.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"

using namespace Tahoe;

/* constructor */
GaoJi::GaoJi(double A, double B, double C, double L_0):
	fA(A),
	fB(B),
	fC(C),
	fL_0(L_0)
{
	SetName("Gao-Ji");
// should insert some consistency checks on the parameters
}

GaoJi::GaoJi(void):
	fA(0.0),
	fB(0.0),
	fC(0.0),
	fL_0(0.0)
{
	SetName("Gao-Ji");
}

/* I/O */
void GaoJi::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      A = " <<   fA << '\n';
	out << "      B = " <<   fB << '\n';
	out << "      C = " <<   fC << '\n';
	out << "    L_0 = " << fL_0 << '\n';
}

void GaoJi::PrintName(ostream& out) const
{
	out << "    Gao-Ji\n";
}

/* returning values */
double GaoJi::Function(double x) const
{
	double dr = x - fL_0;
        double luk;
        double ghos = 1. + dr/(fB*fC);
        double ee = 2.718281828;

        if(ghos < 1./ee)
          luk = -1.;
        else
          luk = log(1. + dr/(fB*fC));

        if(fC > 1.)

	  return( -((fA*fB*fC*(dr - fB*fC*(-1. + pow(1. + dr/(fB*fC),fC)) +
                  (-1. + fC)*(fB*fC + dr)*luk))/
		  (pow(-1. + fC,2.)*pow(1. + dr/(fB*fC),fC))) );
        else
          return 0.0;

}

double GaoJi::DFunction(double x) const
{
	double dr  =  x - fL_0;
        double luk;
        double ghos = 1. + dr/(fB*fC);
        double ee = 2.718281828;

        if(ghos < 1./ee)
          luk = -1.;
        else
          luk = log(1. + dr/(fB*fC));

        if(fC > 1.)
          return( (fA*fB*fC*luk)/pow(1. + dr/(fB*fC),fC) );
        else
          return 0.0;

}

double GaoJi::DDFunction(double x) const
{
	double dr  =  x - fL_0;
        double luk;
        double ghos = 1. + dr/(fB*fC);
        double ee = 2.718281828;

        if(ghos < 1./ee)
          luk = -1.;
        else
          luk = log(1. + dr/(fB*fC));

	if(fC > 1.)
          return( fA*(1. - fC*luk)/pow(1. + dr/(fB*fC),fC+1.) );
        else
          return 0.0;

}

/* returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above */

dArrayT& GaoJi::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
        double luk, ghos;
        double ee = 2.718281828;

	for (int i = 0; i < in.Length(); i++)
	{
		double dr = (*pl++) - fL_0;
                ghos = 1. + dr/(fB*fC);
                
        if(ghos < 1./ee)
          luk = -1.;
        else
          luk = log(1. + dr/(fB*fC));
	
		*pU++ = (fC > 1.) ? -((fA*fB*fC*(dr - fB*fC*(-1. + pow(1. + dr/(fB*fC),fC)) +
                  (-1. + fC)*(fB*fC + dr)*luk))/
				     (pow(-1. + fC,2.)*pow(1. + dr/(fB*fC),fC))) : 0.0;
	}
	return(out);
}

dArrayT& GaoJi::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();
        double luk, ghos;
        double ee = 2.718281828;

	for (int i = 0; i < in.Length(); i++)
	{
		double dr  = (*pl++) - fL_0;
                ghos = 1. + dr/(fB*fC);
		
        if(ghos < 1./ee)
          luk = -1.;
        else
          luk = log(1. + dr/(fB*fC));
	
		*pdU++ = (fC > 1.) ? (fA*fB*fC*luk)/pow(1. + dr/(fB*fC),fC) : 0.0;
	}
	
	return out;
}

dArrayT& GaoJi::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
        double luk, ghos;
        double ee = 2.718281828;

	for (int i = 0; i < in.Length(); i++)
	{
		double dr  = (*pl++) - fL_0;
                ghos = 1. + dr/(fB*fC);
		
        if(ghos < 1./ee)
          luk = -1.;
        else
          luk = log(1. + dr/(fB*fC));
	
		*pddU++ = (fC > 1.) ? fA*(1. - fC*luk)/pow(1. + dr/(fB*fC),fC+1.) : 0.0;
	}
	
	return out;
}

/* describe the parameters needed by the interface */
void GaoJi::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	C1FunctionT::DefineParameters(list);

	// should insert some consistency checks on the parameters
	list.AddParameter(fA, "A");
	list.AddParameter(fB, "B");
	list.AddParameter(fC, "C");
	list.AddParameter(fL_0, "L_0");
}

/* accept parameter list */
void GaoJi::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	C1FunctionT::TakeParameterList(list);

	fA = list.GetParameter("A");
	fB = list.GetParameter("B");
	fC = list.GetParameter("C");
	fL_0 = list.GetParameter("L_0");
}
