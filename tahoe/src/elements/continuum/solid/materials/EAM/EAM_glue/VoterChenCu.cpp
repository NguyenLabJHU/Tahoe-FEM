/* $Id: VoterChenCu.cpp,v 1.5 2004-07-15 08:26:52 paklein Exp $ */
/* created: paklein (12/04/1996) */
#include "VoterChenCu.h"
#include <math.h>
#include "CubicSplineT.h"

using namespace Tahoe;

/* lattice parameters */
const double kLatticeParameterCu = 3.615;  /* angstrom */
const double kCutoffRadiusCu     = 4.961;  /* angstrom */

/* JonnyZ's embedding energy spline data */
#include "VoterChenCu.dat"

/* constructor */
VoterChenCu::VoterChenCu(CBLatticeT& lattice):
	EAM(lattice)
{

}

/* unstressed lattice parameter */
double VoterChenCu::LatticeParameter(void) const
{
	return(kLatticeParameterCu);
}

/**********************************************************************
* Private
**********************************************************************/

void VoterChenCu::SetPairPotential(void)
{	
	fPairPotential = new VCPairPotentialCu();
	if (!fPairPotential) throw ExceptionT::kOutOfMemory;
}

void VoterChenCu::SetElectronDensity(void)
{	
	fElectronDensity = new VCElectronDensityCu();
	if (!fElectronDensity) throw ExceptionT::kOutOfMemory;
}

void VoterChenCu::SetEmbeddingEnergy(void)
{
	/* shallow copy of data */
	dArrayT	  knots(npoints, VC_U_knots_dat);
	dArray2DT coeff(npoints+1, 4, VC_U_coeff_dat);

	fEmbeddingEnergy = new CubicSplineT(knots, coeff);
	if (!fEmbeddingEnergy) throw ExceptionT::kOutOfMemory;
}

/**********************************************************************
*  Glue functions
**********************************************************************/

/* private constructor */
VCPairPotentialCu::VCPairPotentialCu(void)
{
	/* parameters */
	fD		= 0.7366;
	falpha	= 1.919;
	fR		= 2.3250;
}

/* I/O */
void VCPairPotentialCu::Print(ostream& out) const
{
	out << " Parameters:\n\n";
	out << "     D = " << fD     << '\n';
	out << " alpha = " << falpha << '\n';
	out << "     R = " << fR     << '\n';
	out << " r_cut = " << kCutoffRadiusCu  << '\n';

}
	    	   	
void VCPairPotentialCu::PrintName(ostream& out) const
{
	out << "    Voter-Chen Cu pair potential\n";
}     	    	

/* returning values */
double VCPairPotentialCu::Function(double x) const
{
	return( function(x) );	
}

double VCPairPotentialCu::DFunction(double x) const
{
	return( Dfunction(x) );	
}

double VCPairPotentialCu::DDFunction(double x) const
{
	return( DDfunction(x) );	
}

/* returning values in groups */
dArrayT& VCPairPotentialCu::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin = in.Pointer();
	double *pout = out.Pointer();
	int length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = function(*pin++);

	return(out);
}

dArrayT& VCPairPotentialCu::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin = in.Pointer();
	double *pout = out.Pointer();
	int length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = Dfunction(*pin++);

	return(out);
}

dArrayT& VCPairPotentialCu::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin = in.Pointer();
	double *pout = out.Pointer();
	int length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = DDfunction(*pin++);

	return(out);
}

/* return 0th, 1st, and 2nd derivative in the respective
* fields of the dArrayT */  	
void VCPairPotentialCu::SetAll(double x, dArrayT& data) const
{
	data[0] = function(x);
	data[1] = Dfunction(x);
	data[2] = DDfunction(x);
}   	

/**********************************************************************
* Private
**********************************************************************/

/* non-virtual function calls */
double VCPairPotentialCu::function(double x) const
{
	if (x > kCutoffRadiusCu)
		return(0.0);
	else
	{
		double z1, z2, z3, z4, z5;
	
		z1 = -1.*fR;
		z2 = 1.0;
		z3 = pow(x/kCutoffRadiusCu,20);
		z4 = kCutoffRadiusCu + z1;
		z1 = x + z1;
		z2 = -1.*z2*z3;
		z3 = -1.*falpha;
		z2 = 1. + z2;
		z4 = z3*z4;
		z1 = z1*z3;
		z1 = exp(z1);
		z3 = exp(z4);
		z1 = -1.*z1;
		z4 = -1.*z3;
		z1 = 1. + z1;
		z4 = 1. + z4;
		z1 = pow(z1,2.);
		z5 = pow(z4,2.);
		z2 = 0.1*falpha*fD*kCutoffRadiusCu*z2*z3*z4;
		z1 = -1. + z1;
		z3 = -1. + z5;
		z1 = fD*z1;
		z3 = -1.*fD*z3;			

		return z1 + z2 + z3;
	}
}

double VCPairPotentialCu::Dfunction(double x) const
{
	if (x > kCutoffRadiusCu)
		return(0.0);
	else
	{
		double z1, z2, z3, z4, z5, z6, z7;
	
		z1 = -1.*fR;
		z2 = 1.0;
		z3 = pow(x/kCutoffRadiusCu,19.);
		z4 = kCutoffRadiusCu + z1;
		z1 = x + z1;
		z5 = -1.*falpha;
		z4 = z4*z5;
		z1 = z1*z5;
		z4 = exp(z4);
		z1 = exp(z1);
		z5 = -1.*z1;
		z6 = -1.*z4;
		z5 = 1. + z5;
		z6 = 1. + z6;
		z7 = falpha*fD;
		z1 = 2.*z1*z5*z7;
		z2 = -2.*z2*z3*z4*z6*z7;
			
		return z1 + z2;
	}
}

double VCPairPotentialCu::DDfunction(double x) const
{
	if (x > kCutoffRadiusCu)
		return(0.0);
	else
	{
		double z1, z2, z3, z4, z5, z6, z7, z8;

		z1 = pow(falpha,2.);
		z2 = -1.*fR;
		z3 = pow(kCutoffRadiusCu,-1.0);
		z4 = pow(x/kCutoffRadiusCu,18.);
		z5 = kCutoffRadiusCu + z2;
		z2 = x + z2;
		z6 = -1.*falpha;
		z7 = -2.*falpha*z2;
		z5 = z5*z6;
		z2 = z2*z6;
		z2 = exp(z2);
		z6 = exp(z7);
		z5 = exp(z5);
		z7 = -1.*z2;
		z6 = 2.*fD*z1*z6;
		z8 = -1.*z5;
		z7 = 1. + z7;
		z8 = 1. + z8;
		z1 = -2.*fD*z1*z2*z7;
		z2 = -38.*falpha*fD*z3*z4*z5*z8;
	
		return z1 + z2 + z6;
	}
}	

/* private constructor */
VCElectronDensityCu::VCElectronDensityCu(void)
{
	/* parameters */
	fbeta = 4.043;
}

/* I/O */
void VCElectronDensityCu::Print(ostream& out) const
{
	out << " Parameters:\n\n";
	out << "  beta = " << fbeta << '\n';
}
	    	   	
void VCElectronDensityCu::PrintName(ostream& out) const
{
	out << "    Voter-Chen electron density\n";
}     	    	

/* returning values */
double VCElectronDensityCu::Function(double x) const
{
	return( function(x) );	
}

double VCElectronDensityCu::DFunction(double x) const
{
	return( Dfunction(x) );	
}

double VCElectronDensityCu::DDFunction(double x) const
{
	return( DDfunction(x) );	
}

/* returning values in groups */
dArrayT& VCElectronDensityCu::MapFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin = in.Pointer();
	double *pout = out.Pointer();
	int length =  in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = function(*pin++);

	return(out);
}

dArrayT& VCElectronDensityCu::MapDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin = in.Pointer();
	double *pout = out.Pointer();
	int length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = Dfunction(*pin++);

	return(out);
}

dArrayT& VCElectronDensityCu::MapDDFunction(const dArrayT& in, dArrayT& out) const
{
	/* dimension check */
	if ( in.Length() != out.Length() ) throw ExceptionT::kGeneralFail;
	
	const double *pin = in.Pointer();
	double *pout = out.Pointer();
	int length = in.Length();
	
	/* fast mapping */
	for (int i = 0; i < length; i++)
		*pout++ = DDfunction(*pin++);

	return(out);
}

/* return 0th, 1st, and 2nd derivative in the respective
* fields of the dArrayT */  	
void VCElectronDensityCu::SetAll(double x, dArrayT& data) const
{
	data[0] = function(x);
	data[1] = Dfunction(x);
	data[2] = DDfunction(x);
}   	

/**********************************************************************
* Private
**********************************************************************/

/* non-virtual function calls */
double VCElectronDensityCu::function(double x) const
{
	if (x > kCutoffRadiusCu)
		return(0.0);
	else
	{
		double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10;

		z1 = 1.0;
		z2 = -2.*fbeta*kCutoffRadiusCu;
		z3 = -1.*fbeta*kCutoffRadiusCu;
		z4 = pow(kCutoffRadiusCu,5);
		z5 = pow(kCutoffRadiusCu,6);
		z6 = -2.*fbeta*x;
		z7 = -1.*fbeta*x;
		z8 = pow(x,6);
		z9 = pow(x/kCutoffRadiusCu,20);
		z2 = exp(z2);
		z3 = exp(z3);
		z6 = exp(z6);
		z7 = exp(z7);
		z1 = -1.*z1*z9;
		z9 = 512.*z2;
		z2 = -1024.*fbeta*z2;
		z10 = -1.*fbeta*z3;
		z6 = 512.*z6;
		z1 = 1. + z1;
		z3 = z3 + z9;
		z2 = z10 + z2;
		z6 = z6 + z7;
		z4 = 6.*z3*z4;
		z3 = -1.*z3*z5;
		z2 = z2*z5;
		z5 = z6*z8;
		z2 = z2 + z4;
		z1 = 0.05*kCutoffRadiusCu*z1*z2;

		return z1 + z3 + z5;
	}
}

double VCElectronDensityCu::Dfunction(double x) const
{
	if (x > kCutoffRadiusCu)
		return(0.0);
	else
	{
		double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
		double z13, z14;
	
		z1 = 1.0;
		z2 = -2.*fbeta*kCutoffRadiusCu;
		z3 = -1.*fbeta*kCutoffRadiusCu;
		z4 = pow(kCutoffRadiusCu,5.);
		z5 = pow(kCutoffRadiusCu,6.);
		z6 = -2.*fbeta*x;
		z7 = -1.*fbeta*x;
		z8 = pow(x,5.);
		z9 = pow(x,6.);
		z10 = pow(x/kCutoffRadiusCu,19.);
		z2 = exp(z2);
		z3 = exp(z3);
		z6 = exp(z6);
		z7 = exp(z7);
		z11 = 512.*z2;
		z2 = -1024.*fbeta*z2;
		z12 = -1.*fbeta*z3;
		z13 = 512.*z6;
		z6 = -1024.*fbeta*z6;
		z14 = -1.*fbeta*z7;
		z3 = z11 + z3;
		z2 = z12 + z2;
		z7 = z13 + z7;
		z6 = z14 + z6;
		z6 = z6*z9;
		z3 = 6.*z3*z4;
		z2 = z2*z5;
		z4 = 6.*z7*z8;
		z2 = z2 + z3;
		z1 = -1.*z1*z10*z2;

		return z1 + z4 + z6;
	}
}

double VCElectronDensityCu::DDfunction(double x) const
{
	if (x > kCutoffRadiusCu)
		return(0.0);
	else
	{
		double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
		double z13, z14, z15, z16, z17;
	
		z1 = pow(fbeta,2.);
		z2 = pow(kCutoffRadiusCu,-1);
		z3 = -2.*fbeta*kCutoffRadiusCu;
		z4 = -1.*fbeta*kCutoffRadiusCu;
		z5 = pow(kCutoffRadiusCu,5.);
		z6 = pow(kCutoffRadiusCu,6.);
		z7 = -2.*fbeta*x;
		z8 = -1.*fbeta*x;
		z9 = pow(x,4.);
		z10 = pow(x,5.);
		z11 = pow(x,6.);
		z12 = pow(x/kCutoffRadiusCu,18.);
		z3 = exp(z3);
		z4 = exp(z4);
		z7 = exp(z7);
		z8 = exp(z8);
		z13 = 512.*z3;
		z3 = -1024.*fbeta*z3;
		z14 = -1.*fbeta*z4;
		z15 = 512.*z7;
		z16 = -1024.*fbeta*z7;
		z7 = 2048.*z1*z7;
		z17 = -1.*fbeta*z8;
		z1 = z1*z8;
		z4 = z13 + z4;
		z3 = z14 + z3;
		z8 = z15 + z8;
		z13 = z16 + z17;
		z1 = z1 + z7;
		z7 = 30.*z8*z9;
		z8 = 12.*z10*z13;
		z1 = z1*z11;
		z4 = 6.*z4*z5;
		z3 = z3*z6;
		z3 = z3 + z4;
		z2 = -19.*z12*z2*z3;	
		
		return z1 + z2 + z7 + z8;
	}
}	
