/* $Id: ErrorFunc.cpp,v 1.7 2003-11-21 22:41:27 paklein Exp $ */

#include "ErrorFunc.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"
#include "Gamma.h"


using namespace Tahoe;

const int ITMAX = 100;
const double EPS = 3.0e-7;
const double FPMIN = 1.0e-30;

/*
* constructors
*/
ErrorFunc::ErrorFunc()
{
}

/*
* destructors
*/
ErrorFunc::~ErrorFunc()
{
}

/*
* methods
*/
double gammp(double a, double p);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);

/*
* I/O
*/
void ErrorFunc::Print(ostream& out) const
{
	/* parameters */
	out << " No Scaling constant. . . . . . . . . . . . . . . ." << '\n';
}

void ErrorFunc::PrintName(ostream& out) const
{
	out << "    Error Function\n";
}

/*
* Returning values
*/
double ErrorFunc::Function(double x) const
{
	return (x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x));
}

double ErrorFunc::DFunction(double x) const
{
	cout << "\n Derivative of the Error Function not provided!\n";
	throw ExceptionT::kBadInputValue;
	return 0.0*x;	// to avoid generating warning messages
}

double ErrorFunc::DDFunction(double x) const
{
	cout << "\n Second derivative of the Error Function not provided!\n";
	throw ExceptionT::kBadInputValue;
	return 0.0*x;	// to avoid generating warning messages
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& ErrorFunc::MapFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;
	
	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;			
		
		*pU++ = (r < 0.0 ? -gammp(0.5,r*r) : gammp(0.5,r*r));
	}
	return out;
}

dArrayT& ErrorFunc::MapDFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	cout << "\n Derivative of the Bessel Function of the 3rd Kind not tabulated!\n";
	for (int i = 0; i < in.Length(); i++)
	{
//		double r = *pl++;					
	        *pl++;
		*pdU++ = 0.0;
	}
	return out;
}

dArrayT& ErrorFunc::MapDDFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	cout << "\n Second derivative of the Bessel Function of the 3rd Kind not tabulated!\n";
	for (int i = 0; i < in.Length(); i++)
	{
//		double r = *pl++;				
		*pl++;
		*pddU++ = 0.0;
	}
	return out;
}

// Taken from "Numerical Recipes in C", pg. 218-219. This returns the incomplete
// gamma function, P(a,x), evaluated by its series representation.
// Also returns log(Gamma(a)) as gln.
void ErrorFunc::gser(double *gamser, double a, double x, double *gln) const
{
	Gamma gammnln;
	double sum, del, ap;
	
	*gln = log(gammnln.Function(a));
	if (x <= 0.0)
	{
		if (x < 0.0)
		{
			cout << "\n*** x less than 0 in method gser.\n";
			throw ExceptionT::kBadInputValue;
		}
		*gamser = 0.0;
		return;
	}
	else
	{
		ap = a;
		del = sum = 1.0/a;
		for (int n=1; n<=ITMAX; n++)
		{
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS)
			{
				*gamser = sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		cout << "\n*** a too large, ITMAX too small in method gser.\n";
		throw ExceptionT::kBadInputValue;
	}
}

// Taken from "Numerical Recipes in C", pg. 219. This returns the incomplete
// gamma function, Q(a,x), evaluated by its continued fraction representation.
// Also returns log(Gamma(a)) as gln.
void ErrorFunc::gcf(double *gammcf, double a, double x, double *gln) const
{
	Gamma gammnln;
	int i;
	double an, b, c, d, del, h;
	
	*gln = log(gammnln.Function(a));
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	for (i=1; i<=ITMAX; i++)
	{
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = b+an/c;
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (fabs(del-1.0) < EPS)
			break;
	}
	if (i > ITMAX)
	{
		cout << "\n*** a too large, ITMAX too small in method gcf.\n";
		throw ExceptionT::kBadInputValue;
	}
	*gammcf = exp(-x+a*log(x)-(*gln))*h;
}

// Taken from "Numerical Recipes in C", pg. 218. This returns the incomplete
// gamma function, P(a,x).
double ErrorFunc::gammp(double a, double x) const
{
	double gamser, gammcf, gln;
	
	if (x < 0.0 || a <= 0.0)
	{
		cout << "\n*** Invalid arguments in method gammp.\n";
		throw ExceptionT::kBadInputValue;
	}
	if (x < (a+1.0))
	{
		gser(&gamser,a,x,&gln);
		return gamser;
	}
	else
	{
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

