/* $Id: MajumdarBhushan.cpp,v 1.7 2011-12-01 20:37:57 beichuan Exp $ */
#include "MajumdarBhushan.h"
#include <cmath>
#include <iostream>
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

	if (((fD>1.0) && (fD<2.0)) && (fS>0.0) && (fC>0.0))
	{
		ErrorFunc f;
		double hec = 0.5*((2.0-fD)/fD)*(1.0-f.Function(x/(sqrt(2.0)*fS)));
		
		value = pow(fC,1.0-0.5*fD)*pow(hec,0.5*fD)-hec;
		//value = (1.0/(2.0*fD))*(1.0-f.Function(x/(sqrt(2.0)*fS)));	
	}
	else
	{
		if ((fD<=1.0) || (fD>=2.0))	// fD is the fractal dimension of the surface profile
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
	
	if (((fD>1.0) && (fD<2.0)) && (fS>0.0) && (fC>0.0))
	{			
		ErrorFunc f; 
		
		double h = 0.5*(2.0-fD)/fD;
		double ratio = x/(fS*sqrt(2.0));
		double hec = h*(1.0-f.Function(ratio));
			
		if (fD==1.5)
		{
			value = pow(hec,0.75)*log(hec/fC); 
		}
		else
			value = pow(hec,0.5*(3.0-fD))-pow(fC,1.5-fD)*pow(hec,0.5*fD);
	}
	else
	{
		if ((fD<=1.0) || (fD>=2.0))
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

	if (((fD>1.0) && (fD<2.0)) && (fS>0.0) && (fC>0.0))
	{
		ErrorFunc f;
		
		double ratio = x/(fS*sqrt(2.0));
		double h = 0.5*(2.0-fD)/fD;
		double erfc = 1.0-f.Function(ratio);
		double hec = h*erfc;
		double eta = (h/fS)*sqrt(2.0/PI)*exp(-ratio*ratio);
		
		if (fD==1.5)
		{
			double numer = h*exp(-ratio*ratio)*(-4.0+3.0*log(fC/hec));
			double denom = 2.0*fS*sqrt(2.0*PI)*pow(hec,0.25);
			
			value = numer/denom;
		}
		else
		{			
			double c0 = 0.5*fD*pow(fC,1.5-fD)*pow(hec,0.5*(fD-2));
			//double c0 = pow(fC,-fD)*exp(-ratio*ratio)*pow(h,-0.5*fD)*pow(erfc,-1.0-0.5*fD);
			double c1 = 0.5*(3-fD)*pow(hec,0.5*(1-fD));
			//double c1 = (fD-3.0)*pow(fC,fD)*pow(h*erfc,1.5)+fD*pow(fC,1.5)*pow(h*erfc,fD);
		
			value = eta*(c0-c1);	
		}	
	}
	else
	{
		if ((fD<=1.0) || (fD>=2.0))
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

	const double* pl = in.Pointer();
	double* pddU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		r = *pl++;
		
		if (((fD>1.0) && (fD<2.0)) && (fS>0.0) && (fC>0.0))
		{
			ErrorFunc f;
			double hec = 0.5*((2.0-fD)/fD)*(1.0-f.Function(r/(sqrt(2.0)*fS)));
		
			value = pow(fC,1.0-0.5*fD)*pow(hec,0.5*fD)-hec;
			//value = (1.0/(2.0*fD))*(1.0-f.Function(r/(sqrt(2.0)*fS)));		
		}
		else
		{
			if ((fD<=1.0) || (fD>=2.0))
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

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{	
		if (((fD>1.0) && (fD<2.0)) && (fS>0.0) && (fC>0.0))
		{
			r = *pl++;
			ErrorFunc f;
		
			double h = 0.5*(2.0-fD)/fD;
			double ratio = r/(fS*sqrt(2.0));
			double hec = h*(1.0-f.Function(ratio));
			
			if (fD==1.5)
			{
				value = pow(hec,0.75)*log(hec/fC); 
			}
			else
				value = pow(hec,0.5*(3.0-fD))-pow(fC,1.5-fD)*pow(hec,0.5*fD);
		}
		else
		{
			if ((fD<=1.0) || (fD>=2.0))
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

	const double* pl = in.Pointer();
	double* pdU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		if (((fD>1.0) && (fD<2.0)) && (fS>0.0) && (fC>0.0))
		{
			r = *pl++;
			ErrorFunc f;
		
			double ratio = r/(fS*sqrt(2.0));
			double h = 0.5*(2.0-fD)/fD;
			double erfc = 1.0-f.Function(ratio);
			double hec = h*erfc;
			double eta = (h/fS)*sqrt(2.0/PI)*exp(-ratio*ratio);
		
			if (fD==1.5)
			{
				double numer = h*exp(-ratio*ratio)*(-4.0+3.0*log(fC/hec));
				double denom = 2.0*fS*sqrt(2.0*PI)*pow(hec,0.25);
			
				value = numer/denom;
			}
			else
			{
				double c0 = 0.5*fD*pow(fC,1.5-fD)*pow(hec,0.5*(fD-2));
				double c1 = 0.5*(3-fD)*pow(hec,0.5*(1-fD));
		
				value = eta*(c0-c1);	
			}		
		}
		else
		{
			if ((fD<=1.0) || (fD>=2.0))
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

