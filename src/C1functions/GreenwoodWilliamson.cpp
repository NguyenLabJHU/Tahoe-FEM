/* $Id: GreenwoodWilliamson.cpp,v 1.1 2002-01-28 18:41:37 dzeigle Exp $ */

#include "GreenwoodWilliamson.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "ModBessel.h"

/* constants */
const double PI = 2.0*acos(0.0);

/*
* constructors
*/
GreenwoodWilliamson::GreenwoodWilliamson(double MU, double SIGMA, double PARAM):
	fM(MU), fS(SIGMA), fP(PARAM) { }

/*
* I/O
*/
void GreenwoodWilliamson::Print(ostream& out) const
{
	/* parameters */
	out << " Potential parameters:\n";
	out << "      MU = " << fM << '\n';
	out << "      SIGMA = " << fS << '\n';
	out << "      PARAM = " << fP << '\n';
}

void GreenwoodWilliamson::PrintName(ostream& out) const
{
	out << "    Greenwood and Williamson\n";
}

/*
* Returning values
*/
double GreenwoodWilliamson::Function(double x) const
{
	return (0.0*x);
}

double GreenwoodWilliamson::DFunction(double x) const
{
	double f0, f1, f2, xval, diff, yval;
	ModBessel k1(0.25), k3(0.75), k5(1.25);
	
	if (fS==0)
	{
		cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	else
	{
		diff = x - fM;
		xval = pow(diff/(2.0*fS),2.0);
		
		if (diff > 0)
		{
			f0 = fP*sqrt(fS*(x-fM))*exp(-xval);
			f1 = 2.0*(pow(x-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
			f2 = -pow(x-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
	
			yval = f0*(f1 + f2);
		}
		else if (diff == 0)
			yval = 6.09752*fP*pow(fS,3.0); 
		else
		{
			cout << "\n*** Negative approach calculated in GreenwoodWilliamson.cpp\n";
			throw eBadInputValue;
		}
	}
		return yval;
}

double GreenwoodWilliamson::DDFunction(double x) const
{
	double f0, f1, f2, df0, df1, df2, df2sum, xval, diff, yval;
	ModBessel k1(0.25), k3(0.75), k5(1.25), kn1(-0.25), k7(1.75), kn3(-0.75), k9(2.25);
	
	if (fS==0)
	{
		cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	else
	{
		diff = x - fM;
		xval = pow(diff/(2.0*fS),2.0);
		
		if (diff > 0)
		{
			f0 = fP*sqrt(fS*(x-fM))*exp(-xval);
			f1 = 2.0*(pow(x-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
			f2 = -pow(x-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
			
			df0 = -fP*exp(-xval)*(pow(x-fM,2.0)-pow(fS,2.0))/(2.0*pow(fS,1.5)*sqrt(x-fM));
			df1 = (x-fM)*(4.0*k1.Function(xval)-(2.0*xval+1.0)*(kn3.Function(xval)+k5.Function(xval)));
			df2sum = kn1.Function(xval)+k1.Function(xval)+k7.Function(xval)+k9.Function(xval);
			df2 = (x-fM)*(-2.0*(k3.Function(xval)+k5.Function(xval))+xval*df2sum);
	
			yval = df0*(f1+f2)+f0*(df1+df2);
		}
		else if (diff == 0)
/* converges - requires more accurate analysis	*/
			yval = -10.0;
		else
		{
			cout << "\n*** Negative approach calculated in GreenwoodWilliamson.cpp\n";
			throw eBadInputValue;
		}
	}
		return yval;
}

/*
* Returning values in groups - derived classes should define
* their own non-virtual function called within this functon
* which maps in to out w/o requiring a virtual function call
* everytime. Default behavior is just to map the virtual functions
* above.
*/
dArrayT& GreenwoodWilliamson::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl   = in.Pointer();
	double* pddU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		*pddU++ = 0.0;
	}
	return(out);
}

dArrayT& GreenwoodWilliamson::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double f0, f1, f2, xval, diff, yval;
		ModBessel k1(0.25), k3(0.75), k5(1.25);
	
		if (fS==0)
		{
			cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		else
		{
			diff = *pl++ - fM;
			xval = pow(diff/(2.0*fS),2.0);
		
			if (diff > 0)
			{
				f0 = fP*sqrt(fS*(*pl++-fM))*exp(-xval);
				f1 = 2.0*(pow(*pl++-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
				f2 = -pow(*pl++-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
		
				yval = f0*(f1 + f2);
			}
			else if (diff == 0)
				yval = 6.09752*fP*pow(fS,3.0); 
			else
			{
				cout << "\n*** Negative approach calculated in GreenwoodWilliamson.cpp\n";
				throw eBadInputValue;
			}
		}
		
		*pU++ = yval;
	}
	return(out);
}

dArrayT& GreenwoodWilliamson::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	
	for (int i = 0; i < in.Length(); i++)
	{
		double f0, f1, f2, df0, df1, df2, df2sum, xval, diff, yval;
		ModBessel k1(0.25), k3(0.75), k5(1.25), kn1(-0.25), k7(1.75), kn3(-0.75), k9(2.25);
	
		if (fS==0)
		{
			cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		else
		{
			diff = *pl++ - fM;
			xval = pow(diff/(2.0*fS),2.0);
		
			if (diff > 0)
			{
				f0 = fP*sqrt(fS*(*pl++-fM))*exp(-xval);
				f1 = 2.0*(pow(*pl++-fM,2.0) + 2.0*pow(fS,2.0))*k1.Function(xval);
				f2 = -pow(*pl++-fM,2.0)*(k3.Function(xval) + k5.Function(xval));
			
				df0 = -fP*exp(-xval)*(pow(*pl++-fM,2.0)-pow(fS,2.0))/(2.0*pow(fS,1.5)*sqrt(*pl++-fM));
				df1 = (*pl++-fM)*(4.0*k1.Function(xval)-(2.0*xval+1.0)*(kn3.Function(xval)+k5.Function(xval)));
				df2sum = kn1.Function(xval)+k1.Function(xval)+k7.Function(xval)+k9.Function(xval);
				df2 = (*pl++-fM)*(-2.0*(k3.Function(xval)+k5.Function(xval))+xval*df2sum);
	
				yval = df0*(f1+f2)+f0*(df1+df2);
			}
			else if (diff == 0)
/* converges - requires more accurate analysis	*/
				yval = -10.0;
			else
			{
				cout << "\n*** Negative approach calculated in GreenwoodWilliamson.cpp\n";
				throw eBadInputValue;
			}
		}
	
		*pdU++ = yval;
	}
	return(out);
}


