/* $Id: MajumdarBhushan.cpp,v 1.2 2003-05-12 22:01:28 dzeigle Exp $ */
#include "MajumdarBhushan.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionT.h"
#include "dArrayT.h"
#include "ErrorFunc.h"


/* constants */

using namespace Tahoe;

const double PI = 2.0*acos(0.0);

/*
* constructors
*/
MajumdarBhushan::MajumdarBhushan(double FRACDIM, double SIGMA, double C):
fD(FRACDIM), fS(SIGMA), fC(C) { }

/*
* destructors
*/
MajumdarBhushan::~MajumdarBhushan()
{
}

/*
* I/O
*/
void MajumdarBhushan::Print(ostream& out) const
{
	/* parameters */
	out << " MajumdarBhushan parameters:\n";
	out << "      FRACTAL DIMENSION = " << fD << '\n';
	out << "      SIGMA = " << fS << '\n';
	out << "	  C = " << fC << '\n';
}

void MajumdarBhushan::PrintName(ostream& out) const
{
	out << "    Majumdar and Bhushan\n";
}

/*
* Returning values
*/
double MajumdarBhushan::Function(double x) const
// Returns the area value ONLY.
{
	double value=0.0;

	if ((fD>0.0) && (fS>0.0) && (fC>0.0))
	{
		ErrorFunc f;
		value = (1.0/(2.0*fD))*(1.0-f.Function(x/(sqrt(2.0)*fS)));		
	}
	else
	{
		if (fD<=0.0)
		{
			cout << "\n*** Bad DIMENSION value in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
		else if (fS<=0.0)
		{
			cout << "\n*** Bad SIGMA value in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
		else
		{
			cout << "\n*** Error in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
	}
	
	return value;
}

double MajumdarBhushan::DFunction(double x) const
// Returns the load ONLY.
{
	double value=0.0;
	
	if ((fD>0.0) && (fS>0.0) && (fC>0.0))
	{			
		ErrorFunc f;
			
		double c0 = fD/(3.0-2.0*fD);
		double c1 = (2.0-fD)/(2.0*fD);
		double ratio = x/(fS*sqrt(2.0));
			
		double amax = c1*(1.0-f.Function(ratio));
		double amin = fC;
			
		value = pow(amax,fD/2.0)*c0*(pow(amax,1.5-fD)-pow(amin,1.5-fD));
	}
	else
	{
		if (fD<=0.0)
		{
			cout << "\n*** Bad DIMENSION value in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
		else if (fS<=0.0)
		{
			cout << "\n*** Bad SIGMA value in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
		else
		{
			cout << "\n*** Error in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
	}
	
	return value;


}

double MajumdarBhushan::DDFunction(double x) const
// Returns the load gradient ONLY.
{
	double value = 0.0;

	if ((fD>0.0) && (fS>0.0) && (fC>0.0))
	{
		ErrorFunc f;
		
		double ratio = x/(fS*sqrt(2.0)), c0 = fD/(3.0-2.0*fD);
		double amin = fC;
		double amax = ((2.0-fD)/(2.0*fD))*(1.0-f.Function(ratio));
		double damax = ((2.0-fD)/(fS*fD*sqrt(2.0*PI)))*exp(-ratio*ratio);
		
		double v1 = (3.0-fD)*pow(amax,0.5*(1.0-fD));
		double v2 = fD*(pow(amin,1.5-fD)*pow(amax,0.5*fD-1.0));
		
		value = c0*0.5*damax*(v1-v2);		
	}
	else
	{
		if (fD<=0.0)
		{
			cout << "\n*** Bad DIMENSION value in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
		else if (fS<=0.0)
		{
			cout << "\n*** Bad SIGMA value in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
		else
		{
			cout << "\n*** Error in MajumdarBhushan.cpp.\n";
			throw ExceptionT::kBadInputValue;
		}
	}
	
	return value;
}




/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& MajumdarBhushan::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		r = *pl++;
		
		if ((fD>0.0) && (fS>0.0) && (fC>0.0))
		{
			ErrorFunc f;
			value = (1.0/(2.0*fD))*(1.0-f.Function(r/(sqrt(2.0)*fS)));		
		}
		else
		{
			if (fD<=0.0)
			{
				cout << "\n*** Bad DIMENSION value in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
			else if (fS<=0.0)
			{
				cout << "\n*** Bad SIGMA value in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
			else
			{
				cout << "\n*** Error in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
		}
			
		*pddU++ = value;
	}
	return(out);
}

dArrayT& MajumdarBhushan::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{	
		if ((fD>0.0) && (fS>0.0) && (fC>0.0))
		{
			r = *pl++;
			
			ErrorFunc f;
			
			double c0 = fD/(3.0-2.0*fD);
			double c1 = (2.0-fD)/(2.0*fD);
			double ratio = r/(fS*sqrt(2.0));
			
			double amax = c1*(1.0-f.Function(ratio));
			double amin = fC;
			
			value = pow(amax,fD/2.0)*c0*(pow(amax,1.5-fD)-pow(amin,1.5-fD));
		}
		else
		{
			if (fD<=0.0)
			{
				cout << "\n*** Bad DIMENSION value in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
			else if (fS<=0.0)
			{
				cout << "\n*** Bad SIGMA value in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
			else
			{
				cout << "\n*** Error in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
		}
		
		*pU++ = (-value);
	}
	return(out);
}

dArrayT& MajumdarBhushan::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		if ((fD>0.0) && (fS>0.0) && (fC>0.0))
		{
			r = *pl++;
			ErrorFunc f;
		
			double ratio = r/(fS*sqrt(2.0)), c0 = fD/(3.0-2.0*fD);
			double amin = fC;
			double amax = ((2.0-fD)/(2.0*fD))*(1.0-f.Function(ratio));
			double damax = ((2.0-fD)/(fS*fD*sqrt(2.0*PI)))*exp(-ratio*ratio);
		
			double v1 = (3.0-fD)*pow(amax,0.5*(1.0-fD));
			double v2 = fD*(pow(amin,1.5-fD)*pow(amax,0.5*fD-1.0));
		
			value = c0*0.5*damax*(v1-v2);		
		}
		else
		{
			if (fD<=0.0)
			{
				cout << "\n*** Bad DIMENSION value in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
			else if (fS<=0.0)
			{
				cout << "\n*** Bad SIGMA value in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
			else
			{
				cout << "\n*** Error in MajumdarBhushan.cpp.\n";
				throw ExceptionT::kBadInputValue;
			}
		}
		
		*pdU++ = (-value);
	}
	return(out);
}

