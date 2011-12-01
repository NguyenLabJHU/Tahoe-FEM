/* $Id: GWPlastic.cpp,v 1.7 2011-12-01 20:37:57 beichuan Exp $ */
#include "GWPlastic.h"
#include <cmath>
#include <iostream>
#include "ExceptionT.h"
#include "dArrayT.h"
#include "PolyDistributionT.h"


/* A Greenwood-Williamson model for cyclic plastic
 *
 * model parameters:
 * E elastic modulus
 * Y yield value
 * L lengthscale
 * a0 asperity area
 *
 * mu : distribution mean
 * sigma : distribution standard deviation
 *
 * history:
 * dmin : minimum approach
 * Np or Ap : number plastic or area plastic
 */


/* constants */

using namespace Tahoe;

const double PI = 2.0*acos(0.0);
static double BIG = 1.e8;

/*
* constructors
*/
GWPlastic::GWPlastic(double MU, double SIGMA,
double MODULUS, double YIELD, double LENGTHSCALE, double ASPERITYAREA,
double ADHESION_ENERGY, double ADHESION_MODULUS):
		fM(MU), fS(SIGMA),
		fE(MODULUS), fY(YIELD), fL(LENGTHSCALE), fa0(ASPERITYAREA),
		fW(ADHESION_ENERGY),fK(ADHESION_MODULUS),
		fdmin(BIG),fAe(0.0),fAp(0.0) 
{	
		fmoment0 = new PolyDistributionT(0,fM,fS);
		fmoment1 = new PolyDistributionT(1,fM,fS);
		fdc = fL*fY/fE;
		if (fW > 0.0) {
			fda = sqrt(2*fW/fK);
		} else {
			fda = 0.0;
			fK = 0.0;
		}
}

/*
* destructors
*/
GWPlastic::~GWPlastic()
{
		delete fmoment0;
		delete fmoment1;
}

void GWPlastic::ResetParameters(double DMIN)
{
	fdmin=DMIN;
}

/*
* I/O
*/
void GWPlastic::Print(ostream& out) const
{
	/* parameters */
	out << " GWPlastic parameters:\n";
	out << "      MU = " << fM << '\n';
	out << "      SIGMA = " << fS << '\n';
	out << "      MODULUS = " << fM << '\n';
	out << "      YIELD = " << fY << '\n';
	out << "      LENGTHSCALE = " << fL << '\n';
	out << "      ASPERITY AREA = " << fa0 << '\n';
	out << "      ADHESION ENERGY = " << fW << '\n';
	out << "      ADHESION MODULUS = " << fK << '\n';
#if 0
	out << "      ELASTIC AREA = " << fAe << '\n';
	out << "      PLASTIC AREA = " << fAp << '\n';
	out << "      MININUM APPROACH = " << fdmin << '\n';
#endif
}

void GWPlastic::PrintName(ostream& out) const
{
	out << "    Greenwood and Williamson - Plastic \n";
}

double GWPlastic::PlasticArea(double d) const
{
	double value =  fa0*((1-(d+fdc)/fL) *fmoment0->Function(d+fdc)
		  				  		  + 1/fL*fmoment1->Function(d+fdc));
	return value;
}

double GWPlastic::DPlasticArea(double d) const
{
	double value=  fa0*(       (-1/fL)*fmoment0->Function(d+fdc)
				      + (1-(d+fdc)/fL)*fmoment0->DFunction(d+fdc)
		  			            + 1/fL*fmoment1->DFunction(d+fdc));
	return value;
}

/*
* Returning values
*/
double GWPlastic::Function(double d) const
{ // calculate real area of contact
	double Avalue=0.0;
	if (d < fdmin) { // plastic loading
		Avalue = fa0*(fmoment0->Function(d)
				  	 -fmoment0->Function(d+fdc) )
				+PlasticArea(d);
	} else {
		if (d-fdmin < fdc) { // elastic
		   Avalue = fa0*(fmoment0->Function(d)
				  	    -fmoment0->Function(fdmin+fdc) )
				+PlasticArea(fdmin);
		} else { // no contact
			Avalue = 0.0;
		}
	}
	return Avalue;
}

double GWPlastic::DFunction(double d) const
{ // calculate load
	double Fvalue=0.0;
	// calculate plastic area as function of *current* dmin
	if (d < fdmin) { // plastic loading
		Fvalue = fa0*fE/fL*(fmoment1->Function(d)
						-fmoment1->Function(d+fdc) )
				+fY*PlasticArea(d)
				-fa0*fK*(fmoment1->Function(d-fda)
						-fmoment1->Function(d) ) ;
	} else {
		if (d-fdmin < fdc) { // elastic
			Fvalue = fa0*fE/fL*(fmoment1->Function(d)
							-fmoment1->Function(fdmin+fdc) )
				   + fE/fL*(fdc-(d-fdmin))*PlasticArea(fdmin)
				   - fa0*fK*(fmoment1->Function(d-fda)
						    -fmoment1->Function(d) ) ;
		} else { // no contact
			if (d-fdmin < fdc +fda) {
				Fvalue = - fa0*fK*(fmoment1->Function(d-fda)
								  -fmoment1->Function(fdmin+fdc) ) 
						 - fK*(fdc+fda-(d-fdmin))*PlasticArea(fdmin);
			} else {
				Fvalue = 0.0;
			}
		}
	}
	return (-Fvalue);
}

double GWPlastic::DDFunction(double d) const
{ // returns the load gradient
	double dFvalue = 0.0;
	if (d < fdmin) { // plastic loading
		dFvalue = fa0*fE/fL*(fmoment1->DFunction(d)
						- fmoment1->DFunction(d+fdc) )
				+ fY*DPlasticArea(d)
				- fa0*fK*(fmoment1->DFunction(d-fda)
						 -fmoment1->DFunction(d) ) ;
	} else {
		if (d-fdmin < fdc) { // elastic
			dFvalue = fE/fL*fa0*fmoment1->DFunction(d)
					- fE/fL*PlasticArea(fdmin)
				- fa0*fK*(fmoment1->DFunction(d-fda)
						 -fmoment1->DFunction(d) ) ;
		} else { // no contact
			if (d-fdmin < fdc +fda) {
				dFvalue = - fa0*fK*(fmoment1->DFunction(d-fda) )
						 + fK*PlasticArea(fdmin);
			} else {
				dFvalue = 0.0;
			}
		}
	}
	return (-dFvalue);
}




/*
* Returning values in groups - cannot be done with how history 
* variables are currently handled 
*/
dArrayT& GWPlastic::MapFunction(const dArrayT& in, dArrayT& out) const
{
throw ExceptionT::kGeneralFail;
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pU = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		r = *pl++;
		
		value = fS/(sqrt(2.0*PI));
		*pU++ = value;
	}
	return(out);
}

dArrayT& GWPlastic::MapDFunction(const dArrayT& in, dArrayT& out) const
{
throw ExceptionT::kGeneralFail;
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl = in.Pointer();
	double* pF = out.Pointer();
	double r, value = 10.0;
	
	for (int i = 0; i < in.Length(); i++)
	{	
		r = *pl++;
		value = 0.0;
		*pF++ = (-value);
	}
	return(out);
}

dArrayT& GWPlastic::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
throw ExceptionT::kGeneralFail;
	/* dimension checks */
	if (in.Length() != out.Length()) throw ExceptionT::kGeneralFail;

	const double* pl  = in.Pointer();
	double* pdF = out.Pointer();
	double r, value = 0.0;
	
	for (int i = 0; i < in.Length(); i++)
	{
		r = *pl++;
		value = 0.0;
		*pdF++ = (-value);
	}
	return(out);
}

