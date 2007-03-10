/* $Id: InvLangevin.cpp,v 1.2 2007-03-10 01:18:50 tdnguye Exp $ */

#include "InvLangevin.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

InvLangevin::InvLangevin(void) 
{ 
	SetName("inverse_langevin");
}


/*
* Returning values
*/
double InvLangevin::Function(double x) const
{
	if (fabs(x) < 0.84136)
		return(1.31446*tan(1.58986*x)+0.91209*x);
	else if(fabs(x) < 1.0 && fabs(x) > 0.84136)
	{
		int sgn = (floor(x) < 0.0) ? -1 : 1;
		return(1.0/(sgn-x));
	}
	else if (x > 1-kSmall)
		ExceptionT::GeneralFail("InvLangevin::Function", "argument must be less than 1.", x);
}

double InvLangevin::DFunction(double x) const
{
	if (fabs(x) < 0.84136)
	{
		double sec = 1.0/cos(1.58986*x);
		return(0.91209+2.0898073756*sec*sec);
	}
	else if(fabs(x) < 1.0 && fabs(x) > 0.84136)
	{
		int sgn = (floor(x) < 0.0) ? -1 : 1;
		return(1.0/((sgn-x)*(sgn-x)));
	}
	else if (x > 1-kSmall)
		ExceptionT::GeneralFail("InvLangevin::DFunction", "argument must be less than 1.", x);
}


double InvLangevin::DDFunction(double x) const
{
	ExceptionT::GeneralFail("InvLangevin::DDFunction", "Function not defined");

}

/* returning values in groups */
dArrayT& InvLangevin::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pr++;
		if (fabs(x) < 0.84136)
			*pU++ = (1.31446*tan(1.58986*x)+0.91209*x);
		else if(fabs(x) < 1.0 && fabs(x) > 0.84136)
		{
			int sgn = (floor(x) < 0.0) ? -1 : 1;
			*pU++ = (1.0/(sgn-x));
		}
		else if (x > 1-kSmall)
			ExceptionT::GeneralFail("InvLangevin::MapFunction", "argument must be less than 1.", x);
	}
	return(out);
}

dArrayT& InvLangevin::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pr = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double x = *pr++;
		if (fabs(x) < 0.84136)
		{
			double sec = 2.08981/cos(1.58986*x);
			*pdU++ = (0.91209+sec*sec);
		}
		else if(fabs(x) < 1.0 && fabs(x) > 0.84136)
		{
			int sgn = (floor(x) < 0.0) ? -1 : 1;
			*pdU++ = (1.0/((sgn-x)*(sgn-x)));
		}
		else if (x > 1-kSmall)
			ExceptionT::GeneralFail("InvLangevin::MapDFunction", "argument must be less than 1.", x);
	}
	return(out);
}

dArrayT& InvLangevin::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	ExceptionT::GeneralFail("InvLangevin::MapDDFunction", "Function not defined");

}

