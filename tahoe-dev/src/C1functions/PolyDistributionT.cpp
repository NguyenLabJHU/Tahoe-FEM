/* $Id: PolyDistributionT.cpp,v 1.4 2011-12-01 20:37:57 beichuan Exp $ */

#include "PolyDistributionT.h"
#include <iostream>
#include <cmath>
#include "ExceptionT.h"
#include "dArrayT.h"

static const double Pi = 4.0*atan(1.0);

/* constructors */

using namespace Tahoe;

PolyDistributionT::PolyDistributionT(double p, double m, double w): 
fPower(p), fMean(m), fWidth(w)  
{ 
	if (!(fPower==0||fPower==1.0||fPower==1.5)) {
		cout << "\n*** Bad POWER value in PolyDistribution.cpp.\n";
		throw ExceptionT::kBadInputValue;
	}
	if (!(fWidth>0)) {
		cout << "\n*** Bad WIDTH value in PolyDistribution.cpp.\n";
		throw ExceptionT::kBadInputValue;
	}
	if (!(fWidth+fMean>0)) {
		cout << "\n*** WIDTH+MEAN must be positive in PolyDistribution.cpp.\n";
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

/*
* Returning values
*/
double PolyDistributionT::Function(double d) const
{
	double value=0.0;

	if ( d < fMean+fWidth ) {
		if (fPower==0.0) {
			if (d > fMean-fWidth) {
				value = -(pow(d - fMean - fWidth,3)*(3*pow(d,2) - 6*d*fMean + 3*pow(fMean,2) + 9*d*fWidth - 9*fMean*fWidth + 8*pow(fWidth,2)))/(16.*pow(fWidth,5)) ;
			} else { value = 1.0; }
		}
		else if (fPower==1.0) {
			if (d > fMean-fWidth) {
				value = (pow(-d + fMean + fWidth,4)*(pow(d,2) - 2*d*fMean + pow(fMean,2) + 4*d*fWidth - 4*fMean*fWidth + 5*pow(fWidth,2)))/(32.*pow(fWidth,5));
			} else { value = fWidth ;}
		}
		else if (fPower==1.5) {
			if (d > fMean-fWidth) {
				value = (45*pow(d,1.5)*sqrt(Pi)*((8*pow(fMean,4)*pow(-d + fMean + fWidth,2.5))/(15.*pow(d,1.5)*sqrt(Pi)) - (16*pow(fMean,2)*pow(fWidth,2)*pow(-d + fMean + fWidth,2.5))/(15.*pow(d,1.5)*sqrt(Pi)) + (8*pow(fWidth,4)*pow(-d + fMean + fWidth,2.5))/(15.*pow(d,1.5)*sqrt(Pi)) - (32*pow(fMean,3)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (32*fMean*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (16*pow(fMean,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(315.*pow(d,1.5)*sqrt(Pi)) - (16*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(945.*pow(d,1.5)*sqrt(Pi)) - (32*fMean*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(16*pow(d,3) + 40*pow(d,2)*fMean + 70*d*pow(fMean,2) + 105*pow(fMean,3) + 40*pow(d,2)*fWidth + 140*d*fMean*fWidth + 315*pow(fMean,2)*fWidth + 70*d*pow(fWidth,2) + 315*fMean*pow(fWidth,2) + 105*pow(fWidth,3))*sqrt(1 - d/(fMean + fWidth)))/(3465.*pow(d,1.5)*sqrt(Pi)) + (8*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(128*pow(d,4) + 320*pow(d,3)*fMean + 560*pow(d,2)*pow(fMean,2) + 840*d*pow(fMean,3) + 1155*pow(fMean,4) + 320*pow(d,3)*fWidth + 1120*pow(d,2)*fMean*fWidth + 2520*d*pow(fMean,2)*fWidth + 4620*pow(fMean,3)*fWidth + 560*pow(d,2)*pow(fWidth,2) + 2520*d*fMean*pow(fWidth,2) + 6930*pow(fMean,2)*pow(fWidth,2) + 840*d*pow(fWidth,3) + 4620*fMean*pow(fWidth,3) + 1155*pow(fWidth,4))*sqrt(1 - d/(fMean + fWidth)))/(45045.*pow(d,1.5)*sqrt(Pi))))/(64.*pow(fWidth,5));
			} else { 
				value = (1287*pow(fMean,4)*sqrt(fWidth) - 2574*pow(fMean,2)*pow(fWidth,2.5) + 1287*pow(fWidth,4.5) - 1287*pow(fMean,4)*sqrt(fWidth/(fMean + fWidth))*sqrt(fMean + fWidth) - 647*pow(fWidth,4)*sqrt(fWidth/(fMean + fWidth))*sqrt(fMean + fWidth) + 2574*pow(fMean,2)*pow(fWidth/(fMean + fWidth),2.5)*pow(fMean + fWidth,2.5))/(429.*sqrt(2.0)*pow(fWidth,3));
			}
		}
		else {
			cout << "*** ERROR! PolyDistribution p="<<fPower
					<<"  potential unavailable.";
			throw ExceptionT::kBadInputValue;
		}
	}
	
	return value;
}

double PolyDistributionT::DFunction(double d) const
{
	double value=0.0;

	if ( d < fMean+fWidth ) {
		if (fPower==0.0) {
			if (d > fMean-fWidth) {
				value =-(pow(d - fMean - fWidth,3)*(6*d - 6*fMean + 9*fWidth))/(16.*pow(fWidth,5)) - (3*pow(d - fMean - fWidth,2)*(3*pow(d,2) - 6*d*fMean + 3*pow(fMean,2) + 9*d*fWidth - 9*fMean*fWidth + 8*pow(fWidth,2)))/(16.*pow(fWidth,5));
			} else { value = 0.0; }
		}
		else if (fPower==1.0) {
			if (d > fMean-fWidth) {
				value = (pow(-d + fMean + fWidth,4)*(2*d - 2*fMean + 4*fWidth))/(32.*pow(fWidth,5)) - (pow(-d + fMean + fWidth,3)*(pow(d,2) - 2*d*fMean + pow(fMean,2) + 4*d*fWidth - 4*fMean*fWidth + 5*pow(fWidth,2)))/(8.*pow(fWidth,5));
			} else { value = 0.0 ;}
		}
		else if (fPower==1.5) {
			if (d > fMean-fWidth) {


				value = (45*pow(d,1.5)*sqrt(Pi)*((-4*pow(fMean,4)*pow(-d + fMean + fWidth,1.5))/(3.*pow(d,1.5)*sqrt(Pi)) + (8*pow(fMean,2)*pow(fWidth,2)*pow(-d + fMean + fWidth,1.5))/(3.*pow(d,1.5)*sqrt(Pi)) - (4*pow(fWidth,4)*pow(-d + fMean + fWidth,1.5))/(3.*pow(d,1.5)*sqrt(Pi)) - (4*pow(fMean,4)*pow(-d + fMean + fWidth,2.5))/(5.*pow(d,2.5)*sqrt(Pi)) + (8*pow(fMean,2)*pow(fWidth,2)*pow(-d + fMean + fWidth,2.5))/(5.*pow(d,2.5)*sqrt(Pi)) - (4*pow(fWidth,4)*pow(-d + fMean + fWidth,2.5))/(5.*pow(d,2.5)*sqrt(Pi)) + (16*pow(fMean,3)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth))/(105.*pow(d,1.5)*sqrt(Pi)*sqrt(fMean + fWidth)*sqrt(1 - d/(fMean + fWidth))) - (16*fMean*pow(fWidth,2)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth))/(105.*pow(d,1.5)*sqrt(Pi)*sqrt(fMean + fWidth)*sqrt(1 - d/(fMean + fWidth))) - (8*pow(fMean,2)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2)))/(315.*pow(d,1.5)*sqrt(Pi)*sqrt(fMean + fWidth)*sqrt(1 - d/(fMean + fWidth))) + (8*pow(fWidth,2)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2)))/(945.*pow(d,1.5)*sqrt(Pi)*sqrt(fMean + fWidth)*sqrt(1 - d/(fMean + fWidth))) + (16*fMean*pow(-d + fMean + fWidth,2)*(16*pow(d,3) + 40*pow(d,2)*fMean + 70*d*pow(fMean,2) + 105*pow(fMean,3) + 40*pow(d,2)*fWidth + 140*d*fMean*fWidth + 315*pow(fMean,2)*fWidth + 70*d*pow(fWidth,2) + 315*fMean*pow(fWidth,2) + 105*pow(fWidth,3)))/(3465.*pow(d,1.5)*sqrt(Pi)*sqrt(fMean + fWidth)*sqrt(1 - d/(fMean + fWidth))) - (4*pow(-d + fMean + fWidth,2)*(128*pow(d,4) + 320*pow(d,3)*fMean + 560*pow(d,2)*pow(fMean,2) + 840*d*pow(fMean,3) + 1155*pow(fMean,4) + 320*pow(d,3)*fWidth + 1120*pow(d,2)*fMean*fWidth + 2520*d*pow(fMean,2)*fWidth + 4620*pow(fMean,3)*fWidth + 560*pow(d,2)*pow(fWidth,2) + 2520*d*fMean*pow(fWidth,2) + 6930*pow(fMean,2)*pow(fWidth,2) + 840*d*pow(fWidth,3) + 4620*fMean*pow(fWidth,3) + 1155*pow(fWidth,4)))/(45045.*pow(d,1.5)*sqrt(Pi)*sqrt(fMean + fWidth)*sqrt(1 - d/(fMean + fWidth))) - (64*pow(fMean,3)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (64*fMean*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (64*pow(fMean,3)*sqrt(fMean + fWidth)*(-d + fMean + fWidth)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) - (64*fMean*pow(fWidth,2)*sqrt(fMean + fWidth)*(-d + fMean + fWidth)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (16*pow(fMean,3)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(35.*pow(d,2.5)*sqrt(Pi)) - (16*fMean*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(35.*pow(d,2.5)*sqrt(Pi)) + (16*pow(fMean,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(16*d + 20*fMean + 20*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(315.*pow(d,1.5)*sqrt(Pi)) - (16*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(16*d + 20*fMean + 20*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(945.*pow(d,1.5)*sqrt(Pi)) - (32*pow(fMean,2)*sqrt(fMean + fWidth)*(-d + fMean + fWidth)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(315.*pow(d,1.5)*sqrt(Pi)) + (32*pow(fWidth,2)*sqrt(fMean + fWidth)*(-d + fMean + fWidth)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(945.*pow(d,1.5)*sqrt(Pi)) - (8*pow(fMean,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,2.5)*sqrt(Pi)) + (8*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(315.*pow(d,2.5)*sqrt(Pi)) - (32*fMean*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(48*pow(d,2) + 80*d*fMean + 70*pow(fMean,2) + 80*d*fWidth + 140*fMean*fWidth + 70*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(3465.*pow(d,1.5)*sqrt(Pi)) + (64*fMean*sqrt(fMean + fWidth)*(-d + fMean + fWidth)*(16*pow(d,3) + 40*pow(d,2)*fMean + 70*d*pow(fMean,2) + 105*pow(fMean,3) + 40*pow(d,2)*fWidth + 140*d*fMean*fWidth + 315*pow(fMean,2)*fWidth + 70*d*pow(fWidth,2) + 315*fMean*pow(fWidth,2) + 105*pow(fWidth,3))*sqrt(1 - d/(fMean + fWidth)))/(3465.*pow(d,1.5)*sqrt(Pi)) + (16*fMean*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(16*pow(d,3) + 40*pow(d,2)*fMean + 70*d*pow(fMean,2) + 105*pow(fMean,3) + 40*pow(d,2)*fWidth + 140*d*fMean*fWidth + 315*pow(fMean,2)*fWidth + 70*d*pow(fWidth,2) + 315*fMean*pow(fWidth,2) + 105*pow(fWidth,3))*sqrt(1 - d/(fMean + fWidth)))/(1155.*pow(d,2.5)*sqrt(Pi)) + (8*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(512*pow(d,3) + 960*pow(d,2)*fMean + 1120*d*pow(fMean,2) + 840*pow(fMean,3) + 960*pow(d,2)*fWidth + 2240*d*fMean*fWidth + 2520*pow(fMean,2)*fWidth + 1120*d*pow(fWidth,2) + 2520*fMean*pow(fWidth,2) + 840*pow(fWidth,3))*sqrt(1 - d/(fMean + fWidth)))/(45045.*pow(d,1.5)*sqrt(Pi)) - (16*sqrt(fMean + fWidth)*(-d + fMean + fWidth)*(128*pow(d,4) + 320*pow(d,3)*fMean + 560*pow(d,2)*pow(fMean,2) + 840*d*pow(fMean,3) + 1155*pow(fMean,4) + 320*pow(d,3)*fWidth + 1120*pow(d,2)*fMean*fWidth + 2520*d*pow(fMean,2)*fWidth + 4620*pow(fMean,3)*fWidth + 560*pow(d,2)*pow(fWidth,2) + 2520*d*fMean*pow(fWidth,2) + 6930*pow(fMean,2)*pow(fWidth,2) + 840*d*pow(fWidth,3) + 4620*fMean*pow(fWidth,3) + 1155*pow(fWidth,4))*sqrt(1 - d/(fMean + fWidth)))/(45045.*pow(d,1.5)*sqrt(Pi)) - (4*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(128*pow(d,4) + 320*pow(d,3)*fMean + 560*pow(d,2)*pow(fMean,2) + 840*d*pow(fMean,3) + 1155*pow(fMean,4) + 320*pow(d,3)*fWidth + 1120*pow(d,2)*fMean*fWidth + 2520*d*pow(fMean,2)*fWidth + 4620*pow(fMean,3)*fWidth + 560*pow(d,2)*pow(fWidth,2) + 2520*d*fMean*pow(fWidth,2) + 6930*pow(fMean,2)*pow(fWidth,2) + 840*d*pow(fWidth,3) + 4620*fMean*pow(fWidth,3) + 1155*pow(fWidth,4))*sqrt(1 - d/(fMean + fWidth)))/(15015.*pow(d,2.5)*sqrt(Pi))))/(64.*pow(fWidth,5)) + (135*sqrt(d)*sqrt(Pi)*((8*pow(fMean,4)*pow(-d + fMean + fWidth,2.5))/(15.*pow(d,1.5)*sqrt(Pi)) - (16*pow(fMean,2)*pow(fWidth,2)*pow(-d + fMean + fWidth,2.5))/(15.*pow(d,1.5)*sqrt(Pi)) + (8*pow(fWidth,4)*pow(-d + fMean + fWidth,2.5))/(15.*pow(d,1.5)*sqrt(Pi)) - (32*pow(fMean,3)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (32*fMean*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(2*d + 5*fMean + 5*fWidth)*sqrt(1 - d/(fMean + fWidth)))/(105.*pow(d,1.5)*sqrt(Pi)) + (16*pow(fMean,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(315.*pow(d,1.5)*sqrt(Pi)) - (16*pow(fWidth,2)*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(8*pow(d,2) + 20*d*fMean + 35*pow(fMean,2) + 20*d*fWidth + 70*fMean*fWidth + 35*pow(fWidth,2))*sqrt(1 - d/(fMean + fWidth)))/(945.*pow(d,1.5)*sqrt(Pi)) - (32*fMean*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(16*pow(d,3) + 40*pow(d,2)*fMean + 70*d*pow(fMean,2) + 105*pow(fMean,3) + 40*pow(d,2)*fWidth + 140*d*fMean*fWidth + 315*pow(fMean,2)*fWidth + 70*d*pow(fWidth,2) + 315*fMean*pow(fWidth,2) + 105*pow(fWidth,3))*sqrt(1 - d/(fMean + fWidth)))/(3465.*pow(d,1.5)*sqrt(Pi)) + (8*sqrt(fMean + fWidth)*pow(-d + fMean + fWidth,2)*(128*pow(d,4) + 320*pow(d,3)*fMean + 560*pow(d,2)*pow(fMean,2) + 840*d*pow(fMean,3) + 1155*pow(fMean,4) + 320*pow(d,3)*fWidth + 1120*pow(d,2)*fMean*fWidth + 2520*d*pow(fMean,2)*fWidth + 4620*pow(fMean,3)*fWidth + 560*pow(d,2)*pow(fWidth,2) + 2520*d*fMean*pow(fWidth,2) + 6930*pow(fMean,2)*pow(fWidth,2) + 840*d*pow(fWidth,3) + 4620*fMean*pow(fWidth,3) + 1155*pow(fWidth,4))*sqrt(1 - d/(fMean + fWidth)))/(45045.*pow(d,1.5)*sqrt(Pi))))/(128.*pow(fWidth,5));
			} else { value = 0.0 ;}
		}
		else {
			cout << "*** ERROR! PolyDistribution p="<<fPower
					<<"  potential unavailable.";
			throw ExceptionT::kBadInputValue;
		}
	}
	
	return value;
}

double PolyDistributionT::DDFunction(double d) const
{
#pragma unused (d)

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
throw ExceptionT::kGeneralFail;
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pin = in.Pointer();
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
		cout << "*** ERROR! PolyDistribution p="<<fPower<<"  potential unavailable.";
		throw ExceptionT::kBadInputValue;
	}
			
	return(out);
}

dArrayT& PolyDistributionT::MapDFunction(const dArrayT& in, dArrayT& out) const
{
throw ExceptionT::kGeneralFail;
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pin = in.Pointer();
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
throw ExceptionT::kGeneralFail;
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pin = in.Pointer();
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

