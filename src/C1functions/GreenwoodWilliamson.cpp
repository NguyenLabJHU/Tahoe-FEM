/* $Id: GreenwoodWilliamson.cpp,v 1.16 2002-07-02 20:15:09 dzeigle Exp $ */

#include "GreenwoodWilliamson.h"
#include <math.h>
#include <iostream.h>
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "ErrorFunc.h"
#include "Gamma.h"
#include "Kummer.h"
#include "ModBessel.h"


/* constants */

using namespace Tahoe;

const double PI = 2.0*acos(0.0);
const double GWMaxDom = 25.0;

/*
* constructors
*/
GreenwoodWilliamson::GreenwoodWilliamson(double POWER, double MU, double SIGMA):
fP(POWER), fM(MU), fS(SIGMA) { }

/*
* destructors
*/
GreenwoodWilliamson::~GreenwoodWilliamson()
{
}

/*
* I/O
*/
void GreenwoodWilliamson::Print(ostream& out) const
{
	/* parameters */
	out << " GreenwoodWilliamson parameters:\n";
	out << "      MU = " << fM << '\n';
	out << "      SIGMA = " << fS << '\n';
	out << "      POWER = " << fP << '\n';
}

void GreenwoodWilliamson::PrintName(ostream& out) const
{
	out << "    Greenwood and Williamson\n";
}

/*
* Returning values
*/
double GreenwoodWilliamson::Function(double x) const
// Returns the area value (p=1) ONLY.
{
	double value=0.0;

	if (fP==1.0)
	{
		if (fS==0)
		{
			cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		else
		{
			if ((x>fM) || (x<fM))
			{
				ErrorFunc f;
				double expo, errval;
			
				expo = 0.5*(fM-x)*(fM-x)/(fS*fS);
				errval = f.Function((fM-x)/(sqrt(2.0)*fS));
			
				value = 0.5*(fM+exp(-expo)*sqrt(2.0/PI)*fS-x+(fM-x)*errval);
			}
			else if (x==fM)
				value = fS/(sqrt(2.0*PI));
			
		}
	}
	else if (fP==1.5)
	{
		cout << "*** ERROR! Greenwood and Williamson load potential unavailable.";
		throw eBadInputValue;
	}
	else
	{
		cout << "\n*** Bad POWER value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	
	return value;
}

double GreenwoodWilliamson::DFunction(double x) const
// Returns the load (p=3/2) ONLY.
{
	double value=0.0;
	
	if (fP==1.0)
	{
		cout << "*** ERROR! Greenwood and Williamson area gradient unavailable.";
		throw eBadInputValue;
	}
	else if (fP==1.5)
	{
		if (x>fM)
		{
			double term[3], k1, k3, k5, coef, expo;
			ModBessel f1(0.25), f3(0.75), f5(1.25);
		
			for (int i=0; i<3; i++)	// clear entries
				term[i] = 0.0;
			
			coef = -1.0/(8.0*fS*sqrt(PI*(x-fM)));
			expo = 0.25*((fM-x)/fS)*((fM-x)/fS);
			k1 = f1.Function(expo);
			k3 = f3.Function(expo);
			k5 = f5.Function(expo);
		
			term[0] = -(fM-x)*(fM-x)*(k3 + k5);
			term[1] = 2.0*(fM*fM+2.0*fS*fS-2*fM*x+x*x)*k1+term[0];
			term[2] = exp(-expo)*(fM-x)*term[1];
		
			value = coef*term[2];
		}
		else if (x==fM)
		{
			Gamma g;
		
			value = fS*sqrt(fS)*g.Function(5.0/4.0)/(sqrt(PI)*pow(2.0,0.25));
		}
		else if (x<fM)
		{
			double term[4], expo, m5, m7, gn, gp;
			Gamma h;
		
			for (int i=0; i<4; i++)	// clear entries
				term[i] = 0.0;
			
			gn = h.Function(-0.25);
			gp = h.Function(0.25);
			expo = 0.5*((fM-x)/fS)*((fM-x)/fS);
			
			if (expo < GWMaxDom)
			{
				Kummer f5(5.0/4.0,0.5), f7(7.0/4.0,1.5);
				m5 = f5.Function(expo);
				m7 = f7.Function(expo);
			
	
				term[0] = -2.0*sqrt(2.0)*fS*gp*m5;
				term[1] = 3.0*(fM-x)*gn*m7;
				term[2] = exp(-expo)*sqrt(PI*fS)*(term[0]+term[1]);
			}
			else	// domain value too large - will yield error; use anayltic cancellation to evaluate
			{
				double gb5, ga5, gb7, ga7, a5, b5, a7, b7;
				
				a5 = 5.0/4.0;
				b5 = 0.5;
				a7 = 7.0/4.0;
				b7 = 1.5;
				
				gb5 = h.Function(b5);
				ga5 = h.Function(a5);
				gb7 = h.Function(b7);
				ga7 = h.Function(a7);
				
				term[0] = -2.0*sqrt(2.0)*fS*gp*gb5*pow(expo,a5-b5)/ga5;
				term[1] = 3.0*(fM-x)*gn*gb7*pow(expo,a7-b7)/ga7;
				term[2] = sqrt(PI*fS*(fM-x))*(term[0]+term[1]);
			}
			
			term[3] = 2.0*pow(2.0,0.25)*gp*gn;
			
			if (term[3]==0)
			{
				cout << "**Error! - Zero denominator in GreenwoodWilliamson. **\n";
				throw eBadInputValue;
			}
		
			value = term[2]/term[3];
		}
	}
	else
	{
		cout << "\n*** Bad POWER value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	
	return (-value);


}

double GreenwoodWilliamson::DDFunction(double x) const
// Returns the load gradient ONLY.
{
	double value = 0.0;

	if (fP==1.0)
	{
		cout << "*** ERROR! Greenwood and Williamson area gradient unavailable.";
		throw eBadInputValue;
	}
	else if (fP==1.5)
	{
		if (x>fM)
		{
			double term[8], coef, expo, k1, k3, k5, k7, k9, kn1, kn3, sqmux;
			ModBessel f1(0.25), f3(0.75), f5(1.25), f7(1.75), f9(2.25);
			ModBessel fn1(-0.25), fn3(-0.75);
		
			for (int i=0; i<8; i++)	// clear entries
				term[i] = 0.0;
			
			coef = 1.0/(16.0*sqrt(PI*(x-fM))*fS*(x-fM));
			sqmux = fM*fM+2.0*fS*fS-2.0*fM*x+x*x;
			expo = 0.25*((fM-x)/fS)*((fM-x)/fS);
			k1 = f1.Function(expo);
			k3 = f3.Function(expo);
			k5 = f5.Function(expo);
			k7 = f7.Function(expo);
			k9 = f9.Function(expo);
			kn1 = fn1.Function(expo);
			kn3 = fn3.Function(expo);
		
			term[0] = -(fM-x)*(fM-x)*(kn1+k1+k7+k9)/(fS*fS);
			term[1] = -16.0*k1+8.0*(k3+k5);
			term[2] = 2.0*sqmux*(kn3+k5)/(fS*fS);
			term[3] = -0.5*(fM-x)*(fM-x)*(x-fM)*(term[1]+term[2]+term[0]);
			term[4] = -(fM-x)*(fM-x)*(x-fM)*(2.0*sqmux*k1-(fM-x)*(fM-x)*(k3+k5))/(fS*fS);
			term[5] = 2.0*(x-fM)*(2.0*sqmux*k1-(fM-x)*(fM-x)*(k3+k5));
			term[6] = (fM-x)*(2*sqmux*k1-(fM-x)*(fM-x)*(k3+k5));
			term[7] = exp(-expo)*(term[6]+term[5]+term[4]+term[3]);
		
			value = coef*term[7];
		}
		else if (x==fM)
		{
			Gamma g;
			
			value = -3.0*sqrt(fS*PI)/(pow(2.0,5.0/4.0)*g.Function(0.25));
		}
		else if (x<fM)
		{
			double term[4], expo, m5, m7, m9, m11,sqmux, gn, gp;
			Kummer f5(5.0/4.0,0.5), f7(7.0/4.0,1.5), f9(9.0/4.0,1.5), f11(11.0/4.0,2.5);
			Gamma h;
		
			for (int i=0; i<4; i++)	// clear entries
				term[i] = 0.0;
			
			sqmux = (fM-x)*(fM-x)-fS*fS;
			gn = h.Function(-0.25);
			gp = h.Function(0.25);
			expo = 0.5*((fM-x)/fS)*((fM-x)/fS);
			
			if (expo > GWMaxDom)
			{
				double a[4], b[4], ga[4], gb[4], app[4];
				
				a[0] = 5.0/4.0;
				a[1] = 7.0/4.0;
				a[2] = 9.0/4.0;
				a[3] = 11.0/4.0;
				b[0] = 0.5;
				b[1] = 1.5;
				b[2] = 1.5;
				b[3] = 2.5;
				
				for (int i=0; i<4; i++)
				{
					ga[i] = h.Function(a[i]);
					gb[i] = h.Function(b[i]);
					app[i] = gb[i]*pow(expo,a[i]-b[i])/ga[i];
				}

				term[0] = -(fM-x)*(-20.0*fS*gp*app[2]+7.0*sqrt(2.0)*(fM-x)*gn*app[3]);
				term[1] = 8.0*fS*(x-fM)*gp*app[0]+6.0*sqrt(2.0)*sqmux*gn*app[1];
				term[2] = sqrt(PI*(fM-x))*(term[1]+term[0]);				
			}
			else
			{
				m5 = f5.Function(expo);
				m7 = f7.Function(expo);
				m9 = f9.Function(expo);
				m11 = f11.Function(expo);
			
		
				term[0] = -(fM-x)*(-20.0*fS*gp*m9+7.0*sqrt(2.0)*(fM-x)*gn*m11);
				term[1] = 8.0*fS*(x-fM)*gp*m5+6.0*sqrt(2.0)*sqmux*gn*m7;
				term[2] = exp(-expo)*sqrt(PI*(fM-x))*(term[1]+term[0]);
			}
			
			term[3] = pow(2.0,11.0/4.0)*pow(fS,1.5)*sqrt(fM-x)*gn*gp;
			if (term[3]==0)
			{
				cout << "**Error! - Zero denominator in GreenwoodWilliamson. **\n";
				throw eBadInputValue;
			}
		
			value = term[2]/term[3];
		}
	}
	else
	{
		cout << "\n*** Bad POWER value in GreenwoodWilliamson.cpp.\n";
		throw eBadInputValue;
	}
	
	return (-value);
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
	double r, value = -10.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		r = *pl++;
		
		if (fP==1)
		{
			if (fS==0)
			{
				cout << "\n*** Bad SIGMA value in GreenwoodWilliamson.cpp.\n";
				throw eBadInputValue;
			}
			else
			{
				if ((r>fM) || (r<fM))
				{
					ErrorFunc f;
					double expo, errval;
			
					expo = 0.5*(fM-r)*(fM-r)/(fS*fS);
					errval = f.Function((fM-r)/(sqrt(2.0)*fS));
			
					value = 0.5*(fM+exp(-expo)*sqrt(2.0/PI)*fS-r+(fM-r)*errval);
				}
				else if (r==fM)
					value = fS/(sqrt(2.0*PI));
			}
		}
		else if (fP==1.5)
		{
			cout << "*** ERROR! Greenwood and Williamson potential load unavailable.";
			throw eBadInputValue;
		}
		else
		{
			cout << "\n*** Bad POWER value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
			
		*pddU++ = value;
	}
	return(out);
}

dArrayT& GreenwoodWilliamson::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl = in.Pointer();
	double* pU = out.Pointer();
	double r, value=-10.0;
	
	for (int i = 0; i < in.Length(); i++)
	{	
		if (fP==1.0)
		{
			cout << "*** ERROR! Greenwood and Williamson area gradient unavailable.";
			throw eBadInputValue;
		}
		else if (fP==1.5)
		{
			r = *pl++;
		
			if (r>fM)
			{
				double term[3], k1, k3, k5, coef, expo;
				ModBessel f1(0.25), f3(0.75), f5(1.25);
		
				for (int i=0; i<3; i++)	// clear entries
					term[i] = 0.0;
			
				coef = -1.0/(8.0*fS*sqrt(PI*(r-fM)));
				expo = 0.25*((fM-r)/fS)*((fM-r)/fS);
				k1 = f1.Function(expo);
				k3 = f3.Function(expo);
				k5 = f5.Function(expo);
		
				term[0] = -(fM-r)*(fM-r)*(k3 + k5);
				term[1] = 2.0*(fM*fM+2.0*fS*fS-2*fM*r+r*r)*k1+term[0];
				term[2] = exp(-expo)*(fM-r)*term[1];
		
				value = coef*term[2];		
			}
			else if (r==fM)
			{
				Gamma g;
		
				value = fS*sqrt(fS)*g.Function(5.0/4.0)/(sqrt(PI)*pow(2.0,0.25));
			}
			else if (r<fM)
			{
				double term[4], expo, m5, m7, gn, gp;
				Gamma h;
				Kummer f5(5.0/4.0,0.5), f7(7.0/4.0,1.5);
		
				for (int i=0; i<4; i++)	// clear entries
					term[i] = 0.0;
			
				gn = h.Function(-0.25);
				gp = h.Function(0.25);
				expo = 0.5*((fM-r)/fS)*((fM-r)/fS);
			
				if (expo > GWMaxDom)		// geometric cancellation of 
										// Kummer function and exponential function
				{
					double a[2], b[2], ga[2], gb[2];
					double app[2];
				
					a[0] = 5.0/4.0;
					a[1] = 7.0/4.0;
					b[0] = 0.5;
					b[1] = 1.5;
					ga[0] = h.Function(a[0]);
					ga[1] = h.Function(a[1]);
					gb[0] = h.Function(b[0]);
					gb[1] = h.Function(b[1]);
				
					app[0] = gb[0]*pow(expo,a[0]-b[0])/ga[0];
					app[1] = gb[1]*pow(expo,a[1]-b[1])/ga[1];
				
					term[0] = -2.0*sqrt(2.0)*fS*gp*app[0];
					term[1] = 3.0*(fM-r)*gn*app[1];
					term[2] = sqrt(PI*fS*(fM-r))*(term[0]+term[1]);
				}
				else
				{
					m5 = f5.Function(expo);
					m7 = f7.Function(expo);
			
	
					term[0] = -2.0*sqrt(2.0)*fS*gp*m5;
					term[1] = 3.0*(fM-r)*gn*m7;
					term[2] = exp(-expo)*sqrt(PI*fS*(fM-r))*(term[0]+term[1]);
				}
			
				term[3] = 2.0*pow(2.0,0.25)*sqrt(fM-r)*gp*gn;
				if (term[3]==0)
				{
					cout << "**Error! - Zero denominator in GreenwoodWilliamson. **\n";
					throw eBadInputValue;
				}
		
				value = term[2]/term[3];		
			}
		}
		else
		{
			cout << "\n*** Bad POWER value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		
		*pU++ = (-value);
	}
	return(out);
}

dArrayT& GreenwoodWilliamson::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw eGeneralFail;

	double* pl  = in.Pointer();
	double* pdU = out.Pointer();
	double r, value=-1.0e10;
	
	for (int i = 0; i < in.Length(); i++)
	{
		if (fP==1.0)
		{
			cout << "*** ERROR! Greenwood and Williamson area gradient unavailable.";
			throw eBadInputValue;
		}
		else if (fP==1.5)
		{
			r = *pl++;
			
			if (r>fM)
			{
				double term[8], coef, expo, k1, k3, k5, k7, k9, kn1, kn3, sqmux;
				ModBessel f1(0.25), f3(0.75), f5(1.25), f7(1.75), f9(2.25);
				ModBessel fn1(-0.25), fn3(-0.75);
		
				for (int i=0; i<8; i++)	// clear entries
					term[i] = 0.0;
			
				coef = 1.0/(16.0*sqrt(PI*(r-fM))*fS*(r-fM));
				sqmux = fM*fM+2.0*fS*fS-2.0*fM*r+r*r;
				expo = 0.25*((fM-r)/fS)*((fM-r)/fS);
				k1 = f1.Function(expo);
				k3 = f3.Function(expo);
				k5 = f5.Function(expo);
				k7 = f7.Function(expo);
				k9 = f9.Function(expo);
				kn1 = fn1.Function(expo);
				kn3 = fn3.Function(expo);
		
				term[0] = -(fM-r)*(fM-r)*(kn1+k1+k7+k9)/(fS*fS);
				term[1] = -16.0*k1+8.0*(k3+k5);
				term[2] = 2.0*sqmux*(kn3+k5)/(fS*fS);
				term[3] = -0.5*(fM-r)*(fM-r)*(r-fM)*(term[1]+term[2]+term[0]);
				term[4] = -(fM-r)*(fM-r)*(r-fM)*(2.0*sqmux*k1-(fM-r)*(fM-r)*(k3+k5))/(fS*fS);
				term[5] = 2.0*(r-fM)*(2.0*sqmux*k1-(fM-r)*(fM-r)*(k3+k5));
				term[6] = (fM-r)*(2*sqmux*k1-(fM-r)*(fM-r)*(k3+k5));
				term[7] = exp(-expo)*(term[6]+term[5]+term[4]+term[3]);
		
				value = coef*term[7];
			}
			else if (r==fM)
			{
				Gamma g;
			
				value = -3.0*sqrt(fS*PI)/(pow(2.0,5.0/4.0)*g.Function(0.25));
			}
			else if (r<fM)
			{
				double term[4], expo, m5, m7, m9, m11,sqmux, gn, gp;
				Kummer f5(5.0/4.0,0.5), f7(7.0/4.0,1.5), f9(9.0/4.0,1.5), f11(11.0/4.0,2.5);
				Gamma h;
		
				for (int i=0; i<4; i++)	// clear entries
					term[i] = 0.0;
			
				sqmux = (fM-r)*(fM-r)-fS*fS;
				gn = h.Function(-0.25);
				gp = h.Function(0.25);
				expo = 0.5*((fM-r)/fS)*((fM-r)/fS);
			
				if (expo > GWMaxDom)
				{
					double a[4], b[4], ga[4], gb[4], app[4];
				
					a[0] = 5.0/4.0;
					a[1] = 7.0/4.0;
					a[2] = 9.0/4.0;
					a[3] = 11.0/4.0;
					b[0] = 0.5;
					b[1] = 1.5;
					b[2] = 1.5;
					b[3] = 2.5;
				
					for (int i=0; i<4; i++)
					{
						ga[i] = h.Function(a[i]);
						gb[i] = h.Function(b[i]);
						app[i] = gb[i]*pow(expo,a[i]-b[i])/ga[i];
					}
				
					term[0] = -(fM-r)*(-20.0*fS*gp*app[2]+7.0*sqrt(2.0)*(fM-r)*gn*app[3]);
					term[1] = 8.0*fS*(r-fM)*gp*app[0]+6.0*sqrt(2.0)*sqmux*gn*app[1];
					term[2] = sqrt(PI*(fM-r))*(term[1]+term[0]);
				}
				else
				{
					m5 = f5.Function(expo);
					m7 = f7.Function(expo);
					m9 = f9.Function(expo);
					m11 = f11.Function(expo);
			
		
					term[0] = -(fM-r)*(-20.0*fS*gp*m9+7.0*sqrt(2.0)*(fM-r)*gn*m11);
					term[1] = 8.0*fS*(r-fM)*gp*m5+6.0*sqrt(2.0)*sqmux*gn*m7;
					term[2] = exp(-expo)*sqrt(PI*(fM-r))*(term[1]+term[0]);
				}
			
				term[3] = pow(2.0,11.0/4.0)*pow(fS,1.5)*sqrt(fM-r)*gn*gp;
				if (term[3]==0)
				{
					cout << "**Error! - Zero denominator in GreenwoodWilliamson. **\n";
					throw eBadInputValue;
				}
		
				value = term[2]/term[3];
			}
		}
		else
		{
			cout << "\n*** Bad POWER value in GreenwoodWilliamson.cpp.\n";
			throw eBadInputValue;
		}
		
		*pdU++ = (-value);
	}
	return(out);
}

