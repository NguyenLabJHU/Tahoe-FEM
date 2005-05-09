/* $Id: FungwRep.cpp,v 1.2 2005-05-09 17:19:14 thao Exp $ */

#include "FungwRep.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

FungwRep::FungwRep(double A, double B, double C, double n): 
	fA(A), 
	fB(B),
	fC(C),
	fn(n) 
{ 
	SetName("Fung-Type");
}

/* I/O */
void FungwRep::Print(ostream& out) const
{
	/* parameters */
	out << " fA. . . . . . . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " fB. . . . . . . . . . . . . . . . . . . . . . = " << fB << '\n';
}

void FungwRep::PrintName(ostream& out) const
{
	out << "    Worm Like Chain Statistics\n";
}

/*
* Returning values
*/
double FungwRep::Function(double r) const
{
	return ( fA*(exp(fB*(r*r-1.0))-1.0)+ fC*pow(r, -fn) );
}

double FungwRep::DFunction(double r) const
{
	return ( 2.0*fA*fB*r*exp(fB*(r*r-1.0))- fC*fn*pow(r,-(fn+1.0)) );
}

double FungwRep::DDFunction(double r) const
{
	return ( 2.0*fA*fB*(1.0+2.0*fB*r*r)*exp(fB*(r*r-1.0))+ fC*fn*(fn+1.0)*pow(r, -(fn+2.0)) ) ;
}

/* returning values in groups */
dArrayT& FungwRep::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pU++ = (fA*(exp(fB*(r*r-1.0))-1.0)+ fC*pow(r, -fn));
	}
	return(out);
}

dArrayT& FungwRep::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pdU++ = ( 2.0*fA*fB*r*exp(fB*(r*r-1.0))- fC*fn*pow(r,-(fn+1.0)) );
	}
	return(out);
}

dArrayT& FungwRep::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		*pddU++ = ( 2.0*fA*fB*(1.0+2.0*fB*r*r)*exp(fB*(r*r-1.0))+ fC*fn*(fn+1.0)*pow(r, -(fn+2.0)) );
	}
	return(out);
}

