/* $Id: PolyDistributionT.cpp,v 1.1 2003-06-03 16:32:12 rjones Exp $ */

#include "PolyDistributionT.h"
#include <iostream.h>
#include <math.h>
#include "ExceptionT.h"
#include "dArrayT.h"

/* constructors */

using namespace Tahoe;

PolyDistributionT::PolyDistributionT(double p, double m, double w): 
fPower(p), fMean(m), fWidth(w)  
{ 
	if (!(fPower==0||fPower==1.0||fPower==1.5)) {
		cout << "\n*** Bad POWER value in PolyDistribution.cpp.\n";
		throw ExceptionT::kBadInputValue;
	}
	if (fWidth==0) {
		cout << "\n*** Bad WIDTH value in PolyDistribution.cpp.\n";
		throw ExceptionT::kBadInputValue;
	}
}

/* I/O */
void PolyDistributionT::Print(ostream& out) const
{
	/* parameters */
	out << " Moment. . . . . . . . . . . . . . . . . . . . . = " << fPower << '\n';
	out << " Mean. . . . . . . . . . . . . . . . . . . . . . = " << fMean << '\n';
	out << " Width . . . . . . . . . . . . . . . . . . . . . = " << fWidth << '\n';
}

void PolyDistributionT::PrintName(ostream& out) const
{
	out << "    Polynominal distribution function\n";
}

double PolyDistributionT::Mom0(const double x, const double d) const
{
		double value =(
		((15*pow(fMean,4) - 30*pow(fMean,2)*pow(fWidth,2) + 15*pow(fWidth,4))*x)
		+ ((-30*pow(fMean,3) + 30*fMean*pow(fWidth,2))*pow(x,2))
		+ ((30*pow(fMean,2) - 10*pow(fWidth,2))*pow(x,3))
		+ (3*pow(x,5))
		) /(16.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::dMom0dx(const double x, const double d) const
{
		double value = (
		(15*(pow(fMean,4) - 2*pow(fMean,2)*pow(fWidth,2) + pow(fWidth,4)))
		+ (15*(-4*pow(fMean,3) + 4*fMean*pow(fWidth,2))*x)
		+ (15*(6*pow(fMean,2) - 2*pow(fWidth,2))*pow(x,2))
		- (60*fMean*pow(x,3))
		+ (15*pow(x,4))
		)/(16.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::Mom1(const double x, const double d) const
{
		double value =(
		((-30*d*pow(fMean,4) + 60*d*pow(fMean,2)*pow(fWidth,2) 
			- 30*d*pow(fWidth,4))*x)
		+ ((60*d*pow(fMean,3) + 15*pow(fMean,4) - 60*d*fMean*pow(fWidth,2) 
			- 30*pow(fMean,2)*pow(fWidth,2) + 15*pow(fWidth,4))*pow(x,2))
		+ ((-60*d*pow(fMean,2) - 40*pow(fMean,3) + 20*d*pow(fWidth,2) 
			+ 40*fMean*pow(fWidth,2))*pow(x,3))
		+ ((30*d*fMean + 45*pow(fMean,2) - 15*pow(fWidth,2))*pow(x,4))
		+ ((-6*d - 24*fMean)*pow(x,5))
		+ (5*pow(x,6))
		)/(32.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::dMom1dx(const double x, const double d) const
{
		double value =(
		(-15*(d*pow(fMean,4) - 2*d*pow(fMean,2)*pow(fWidth,2) + d*pow(fWidth,4)))
		- (15*(-4*d*pow(fMean,3) - pow(fMean,4) + 4*d*fMean*pow(fWidth,2) 
		+ 2*pow(fMean,2)*pow(fWidth,2) - pow(fWidth,4))*x)
		- (15*(6*d*pow(fMean,2) + 4*pow(fMean,3) - 2*d*pow(fWidth,2) 
		- 4*fMean*pow(fWidth,2))*pow(x,2))
		- (15*(-4*d*fMean - 6*pow(fMean,2) + 2*pow(fWidth,2))*pow(x,3))
		- (15*(d + 4*fMean)*pow(x,4))
		+ (15*pow(x,5))
		)/(16.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::dMom1dd(const double x, const double d) const
{
		double value =(
		(-15*d*pow(pow(fMean,2) - pow(fWidth,2) - 2*fMean*x + pow(x,2),2))
		+ (15*x*pow(pow(fMean,2) - pow(fWidth,2) - 2*fMean*x + pow(x,2),2))
		)/(16.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::Mom1_5(const double x, const double d) const
{
		double value =(
		((384*pow(d,4) - 2496*pow(d,3)*fMean + 6864*pow(d,2)*pow(fMean,2) 
		- 10296*d*pow(fMean,3) + 9009*pow(fMean,4) - 2288*pow(d,2)*pow(fWidth,2) 
		+ 10296*d*fMean*pow(fWidth,2) - 18018*pow(fMean,2)*pow(fWidth,2) 
		+ 9009*pow(fWidth,4))*pow(-d + x,2.5))
		+ ((960*pow(d,3) - 6240*pow(d,2)*fMean + 17160*d*pow(fMean,2) 
		- 25740*pow(fMean,3) - 5720*d*pow(fWidth,2) 
		+ 25740*fMean*pow(fWidth,2))*x*pow(-d + x,2.5))
		+ ((1680*pow(d,2) - 10920*d*fMean + 30030*pow(fMean,2) 
		- 10010*pow(fWidth,2))*pow(x,2)*pow(-d + x,2.5))
		+ ((2520*d - 16380*fMean)*pow(x,3)*pow(-d + x,2.5))
		)/(24024.*pow(fWidth,5)) 
		+ (15*pow(x,4)*pow(-d + x,2.5))/(104.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::dMom1_5dx(const double x, const double d) const
{
		double value =(
		(15*(pow(fMean,4) - 2*pow(fMean,2)*pow(fWidth,2) + pow(fWidth,4))*pow(-d + x,1.5))
		+ (15*(-4*pow(fMean,3) + 4*fMean*pow(fWidth,2))*x*pow(-d + x,1.5))
		+ (15*(6*pow(fMean,2) - 2*pow(fWidth,2))*pow(x,2)*pow(-d + x,1.5))
		- 4.*(15*fMean*pow(x,3)*pow(-d + x,1.5))
		+ (15*pow(x,4)*pow(-d + x,1.5))
		)/(16.*pow(fWidth,5));
		return value;
}

double PolyDistributionT::dMom1_5dd(const double x, const double d) const
{
		double value =(
		(15*pow(-d + x,1.5)*pow(pow(fMean,2) - pow(fWidth,2) - 2*fMean*x + pow(x,2),2))
		)/(16.*pow(fWidth,5));
		return value;
}

/*
* Returning values
*/
double PolyDistributionT::Function(double d) const
// Returns the area value (p = 1) ONLY.
{
	double value=0.0;

	if (fPower==0.0) {
		double uplimit= (d > fMean+fWidth) ? 0.0      : Mom0(fMean+fWidth,d);
		double lolimit= (d > fMean-fWidth) ? Mom0(d,d): Mom0(fMean-fWidth,d); 
		value = uplimit - lolimit;
	}
	else if (fPower==1.0) {
		double uplimit= (d > fMean+fWidth) ? 0.0      : Mom1(fMean+fWidth,d);
		double lolimit= (d > fMean-fWidth) ? Mom1(d,d): Mom1(fMean-fWidth,d); 
		value = uplimit - lolimit;
	}
	else if (fPower==1.5) {
		double uplimit= (d > fMean+fWidth) ? 0.0        : Mom1_5(fMean+fWidth,d);
		double lolimit= (d > fMean-fWidth) ? Mom1_5(d,d): Mom1_5(fMean-fWidth,d); 
		value = uplimit - lolimit;
	}
	else {
		cout << "*** ERROR! PolyDistribution p="<<fPower<<"  potential unavailable.";
		throw ExceptionT::kBadInputValue;
	}
	
	return value;
}

double PolyDistributionT::DFunction(double d) const
{
	double value=0.0;
	
	if (fPower==0.0) {
	}
	else if (fPower==1.0) {
	}
	if (fPower==1.5) {
	}
	else {
		cout << "*** ERROR! PolyDistribution p="<<fPower<<"  gradient unavailable.";
		throw ExceptionT::kBadInputValue;
	}
	
	return value;


}

double PolyDistributionT::DDFunction(double d) const
{
	double value = 0.0;

	if (0) {
	}
	else
	{
		cout << "*** ERROR! PolyDistribution p="<<fPower<<" 2nd  gradient unavailable.";
		throw ExceptionT::kBadInputValue;
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
dArrayT& PolyDistributionT::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pin  = in.Pointer();
	double* pout = out.Pointer();
	double x, value = 0.0;
	if (fPower==0.0) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			double uplimit= (x > fMean+fWidth) ? 0.0      : Mom0(fMean+fWidth,x);
			double lolimit= (x > fMean-fWidth) ? Mom0(x,x): Mom0(fMean-fWidth,x); 
			value = uplimit - lolimit;
			*pout++ = value;
		}
	}
	else if (fPower==1.0) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			double uplimit= (x > fMean+fWidth) ? 0.0      : Mom1(fMean+fWidth,x);
			double lolimit= (x > fMean-fWidth) ? Mom1(x,x): Mom1(fMean-fWidth,x); 
			value = uplimit - lolimit;
		}
	}
	else if (fPower==1.5) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			double uplimit= (x > fMean+fWidth) ? 0.0        : Mom1_5(fMean+fWidth,x);
			double lolimit= (x > fMean-fWidth) ? Mom1_5(x,x): Mom1_5(fMean-fWidth,x); 
			value = uplimit - lolimit;
			*pout++ = value;
		}
	}
	else {
		cout << "*** ERROR! PolyDistribution p="<<fPower<<"  potential unavailable.";
		throw ExceptionT::kBadInputValue;
	}
			
	return(out);
}

dArrayT& PolyDistributionT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pin  = in.Pointer();
	double* pout = out.Pointer();
	double x, value = 0.0;
	if (fPower==0.0) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			*pout++ = value;
		}
	}
	else if (fPower==1.0) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			*pout++ = value;
		}
	}
	else if (fPower==1.5) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			*pout++ = value;
		}
	}
	else {
		cout << "*** ERROR! PolyDistribution p="<<fPower<<"  gradient unavailable.";
		throw ExceptionT::kBadInputValue;
	}
			
	return(out);
}

dArrayT& PolyDistributionT::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	double* pin  = in.Pointer();
	double* pout = out.Pointer();
	double x, value = 0.0;
	if (0) {
		for (int i = 0; i < in.Length(); i++) {
			x = *pin++;
			*pout++ = value;
		}
	}
	else {
		cout << "*** ERROR! PolyDistribution p="<<fPower<<" 2nd gradient unavailable.";
		throw ExceptionT::kBadInputValue;
	}
			
	return(out);
}

