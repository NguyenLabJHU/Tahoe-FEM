/* $Id: WLCwRep.cpp,v 1.1 2005-04-20 23:47:27 thao Exp $ */

#include "WLCwRep.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

WLCwRep::WLCwRep(double A, double B, double C, double n): 
	fA(A), 
	fB(B),
	fC(C),
	fn(n) 
{ 
	SetName("Smith-Ferrante");
}

/* I/O */
void WLCwRep::Print(ostream& out) const
{
	/* parameters */
	out << " fA. . . . . . . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " fB. . . . . . . . . . . . . . . . . . . . . . = " << fB << '\n';
}

void WLCwRep::PrintName(ostream& out) const
{
	out << "    Worm Like Chain Statistics\n";
}

/*
* Returning values
*/
double WLCwRep::Function(double r) const
{
	double p = r/fB;
	return (fA * (0.5*r*p - 0.25*r + 0.25*fB/(1-p)) + fC*pow(r, -fn));
}

double WLCwRep::DFunction(double r) const
{
	double p = r/fB;
	return (fA * ( p - 0.25 + 0.25/((1-p)*(1-p)) ) - fC*fn*pow(r,-(fn+1.0)) );
}

double WLCwRep::DDFunction(double r) const
{
	double p = r/fB;
	return (fA/fB * ( 1 + 0.5/((1-p)*(1-p)*(1-p)) ) + fC*fn*(fn+1.0)*pow(r, -(fn+2.0)) );
}

/* returning values in groups */
dArrayT& WLCwRep::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		double p = r/fB;
		*pU++ = (fA * (0.5*r*p - 0.25*r + 0.25*fB/(1-p)) + fC*pow(r, -fn) );
	}
	return(out);
}

dArrayT& WLCwRep::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		double p = r/fB;
		*pdU++ = (fA * ( p - 0.25 + 0.25/((1-p)*(1-p)) ) - fC*fn*pow(r,-(fn+1.0)) );
	}
	return(out);
}

dArrayT& WLCwRep::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		double p = r/fB;
		*pddU++ = (fA/fB * ( 1 + 0.5/((1-p)*(1-p)*(1-p)) ) + fC*fn*(fn+1.0)*pow(r, -(fn+2.0)) );
	}
	return(out);
}

