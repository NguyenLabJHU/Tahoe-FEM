/* $Id: Gamma.cpp,v 1.1 2002-04-22 16:09:23 dzeigle Exp $ */
/* created: dzeigle (4/18/2002) */

#include "Gamma.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"

/*
* Constructor
*/
Gamma::Gamma()
{
}

/*
* Destructor
*/
Gamma::~Gamma()
{
}
	 
/*
* I/O
*/
void Gamma::Print(ostream& out) const
{
	/* parameters */
	out << " Gamma function. . . . . . . . . . . . . . . . \n";
}

void Gamma::PrintName(ostream& out) const
{
	out << "    Modified Bessel Function\n";
}

/*
* Returning values
*/
double Gamma::Function(double xx) const
{
	double x, y, tmp, ser;
	static double cof[6] = {76.18009172947146, -86.50532032941677,
	24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for (int j=0; j<=5; j++)
	{
		ser += cof[j]/++y;
	}
	return exp(-tmp+log(2.5066282746310005*ser/x));
}

double Gamma::DFunction(double x) const
{
	cout << "\n First derivative of the Gamma Function not tabulated!\n";
	return 0*x;
}

double Gamma::DDFunction(double x) const
{
	cout << "\n Second derivative of the Gamma Function not tabulated!\n";
	return 0.0*x;
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& Gamma::MapFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;
	
	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double x, y, tmp, ser;
	static double cof[6] = {76.18009172947146, -86.50532032941677,
	24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
	
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;
		
		y = x = r;
		tmp = x+5.5;
		tmp -= (x+0.5)*log(tmp);
		ser = 1.000000000190015;
		for (int j=0; j<=5; j++)
		{
			ser += cof[j]/++y;
		}
		*pU++ =  exp(-tmp+log(2.5066282746310005*ser/x));
	}
	return out;
}

dArrayT& Gamma::MapDFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	cout << "\n Derivative of the Gamma Function not tabulated!\n";
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;					
		*pdU++ = 0.0;
	}
	return out;
}

dArrayT& Gamma::MapDDFunction(const dArrayT& in,  dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	cout << "\n Second derivative of the Gamma Function not tabulated!\n";
	for (int i = 0; i < in.Length(); i++)
	{
		double r = *pl++;				
		*pddU++ = 0.0;
	}
	return out;
}



