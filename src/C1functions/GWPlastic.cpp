/* $Id: GWPlastic.cpp,v 1.2 2003-06-30 22:07:25 rjones Exp $ */
#include "GWPlastic.h"
#include <math.h>
#include <iostream.h>
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
 * Np or Ap : number plastic or area plastic
 * gdot : loading direction
 */


/* constants */

using namespace Tahoe;

const double PI = 2.0*acos(0.0);

/*
* constructors
*/
GWPlastic::GWPlastic(double MU, double SIGMA,
double MODULUS, double YIELD, double LENGTHSCALE, double ASPERITYAREA):
		fM(MU), fS(SIGMA),
		fE(MODULUS), fY(YIELD), fL(LENGTHSCALE), fa0(ASPERITYAREA),
		fdmin(1.e8),fddot(-1),fAe(0.0),fAp(0.0) 
{	
		fmoment0 = new PolyDistributionT(0,fM,fS);
		fmoment1 = new PolyDistributionT(1,fM,fS);
		fdc = fY/fE;
}

/*
* destructors
*/
GWPlastic::~GWPlastic()
{
		delete fmoment0;
		delete fmoment1;
}

void GWPlastic::ResetParameters(double DMIN, double DDOT)
{
	fdmin=DMIN;
	fddot=DDOT;	
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
#if 0
	out << "      GDOT = " << fddot << '\n';
	out << "      ELASTIC AREA = " << fAe << '\n';
	out << "      PLASTIC AREA = " << fAp << '\n';
	out << "      MININUM APPROACH = " << fdmin << '\n';
#endif
}

void GWPlastic::PrintName(ostream& out) const
{
	out << "    Greenwood and Williamson - Plastic \n";
}

/*
* Returning values
*/
double GWPlastic::Function(double d) const
{ // calculate real area of contact
	double value=0.0;
	return value;
}

double GWPlastic::DFunction(double d) const
{ // calculate load
	double Fvalue=0.0;
	double Ap = fa0*((1-(fdmin+fdc)/fL) *fmoment0->Function(fdmin+fdc)
					+ 1/fL*fmoment1->Function(fdmin+fdc));
	// calculate plastic area as function of *current* dmin
	if (fddot < 0.0 && d < fdmin) { // plastic loading
		Fvalue = fa0*fE*(fmoment1->Function(d)
						-fmoment1->Function(d+fdc) )
				+fY*Ap;
	} else {
		if (d-fdmin < fdc) { // elastic
			Fvalue = fE*fa0*(fmoment1->Function(d)
							-fmoment1->Function(fdmin+fdc) )
					+fE*(fdc-(d-fdmin))*Ap;
		} else { // no contact
			Fvalue = 0.0;
		}
	}
	return (-Fvalue);
}

double GWPlastic::DDFunction(double d) const
{ // returns the load gradient
	double dFvalue = 0.0;
	double Ap = fa0*((1-(fdmin+fdc)/fL) *fmoment0->Function(fdmin+fdc)
					+ 1/fL*fmoment1->Function(fdmin+fdc));
	if (fddot < 0.0 && d < fdmin) { // plastic loading
		dFvalue = fa0*fE*(fmoment1->DFunction(d)
						- fmoment1->DFunction(d+fdc) )
				+fY*Ap;
	} else {
		if (d-fdmin < fdc) { // elastic
			dFvalue = fE*fa0*(fmoment1->DFunction(d)
							- fmoment1->DFunction(fdmin+fdc))
					-fE*d*Ap;
		} else { // no contact
			dFvalue = 0.0;
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

	double* pl = in.Pointer();
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

	double* pl = in.Pointer();
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

	double* pl  = in.Pointer();
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

