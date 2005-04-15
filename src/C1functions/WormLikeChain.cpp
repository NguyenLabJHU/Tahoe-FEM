/* $Id: WormLikeChain.cpp,v 1.1 2005-04-15 22:47:45 thao Exp $ */

#include "WormLikeChain.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

WormLikeChain::WormLikeChain(double A, double B): 
	fA(A), 
	fB(B) 
{ 
	SetName("Smith-Ferrante");
}

/* I/O */
void WormLikeChain::Print(ostream& out) const
{
	/* parameters */
	out << " fA. . . . . . . . . . . . . . . . . . . . . = " << fA << '\n';
	out << " fB. . . . . . . . . . . . . . . . . . . . . . = " << fB << '\n';
}

void WormLikeChain::PrintName(ostream& out) const
{
	out << "    Worm Like Chain Statistics\n";
}

/*
* Returning values
*/
double WormLikeChain::Function(double r) const
{
	double p = r/fB;
	return (fA * (0.5*r*p - 0.25*r + 0.25*fB/(1-p)) );
}

double WormLikeChain::DFunction(double r) const
{
	double p = r/fB;
	return (fA * ( p - 0.25 + 0.25/((1-p)*(1-p)) ) );
}

double WormLikeChain::DDFunction(double r) const
{
	double p = r/fB;
	return (fA/fB * ( 1 + 0.5/((1-p)*(1-p)*(1-p)) ) );
}

/* returning values in groups */
dArrayT& WormLikeChain::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		double p = r/fB;
		*pU++ = (fA * (0.5*r*p - 0.25*r + 0.25*fB/(1-p)));
	}
	return(out);
}

dArrayT& WormLikeChain::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		double p = r/fB;
		*pdU++ = (fA * ( p - 0.25 + 0.25/((1-p)*(1-p)) ) );
	}
	return(out);
}

dArrayT& WormLikeChain::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = (*pr++);
		double p = r/fB;
		*pddU++ = (fA/fB * ( 1 + 0.5/((1-p)*(1-p)*(1-p)) ) );
	}
	return(out);
}

