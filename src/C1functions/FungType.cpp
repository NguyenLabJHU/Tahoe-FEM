/* $Id: FungType.cpp,v 1.1 2005-04-15 22:47:45 thao Exp $ */

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
	SetName("Fung-Type");
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
	out << "    Worm Like Chain Statistics\n";
}

/*
* Returning values
*/
double FungType::Function(double r) const
{
	return ( fA*(exp(fB*(r-1.0))-1.0) );
}

double FungType::DFunction(double r) const
{
	return ( fA*fB*exp(fB*(r-1.0)) );
}

double FungType::DDFunction(double r) const
{
	return ( fA*fB*fB*exp(fB*(r-1.0)) );
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
		*pU++ = (fA*(exp(fB*(r-1.0))-1.0));
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
		*pdU++ = ( fA*fB*exp(fB*(r-1.0)) );
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
		*pddU++ = (fA*fB*fB*exp(fB*(r-1.0)) );
	}
	return(out);
}

