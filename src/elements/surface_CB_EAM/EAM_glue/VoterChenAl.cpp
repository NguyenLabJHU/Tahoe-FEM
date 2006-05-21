/* $Id: VoterChenAl.cpp,v 1.1 2006-05-21 15:55:19 hspark Exp $ */
/* created: paklein (12/04/1996) */
#include "VoterChenAl.h"
#include <math.h>
#include "CubicSplineT.h"

using namespace Tahoe;

/* lattice parameters */
const double kLatticeParameterAl = 4.05; /* angstrom */
const double kCutoffRadiusAl 	 = 5.5550; /* angstrom */

/* constructor */
VoterChenAl::VoterChenAl(CBLatticeT& lattice):
	EAM(lattice)
{

}

/* unstressed lattice parameter */
double VoterChenAl::LatticeParameter(void) const
{
	return(kLatticeParameterAl);
}

/**********************************************************************
* Private
**********************************************************************/

void VoterChenAl::SetPairPotential(void)
{	
	fPairPotential = new VCPairPotential();
	if (!fPairPotential) throw ExceptionT::kOutOfMemory;
}

void VoterChenAl::SetElectronDensity(void)
{	
	fElectronDensity = new VCElectronDensity();
	if (!fElectronDensity) throw ExceptionT::kOutOfMemory;
}

void VoterChenAl::SetEmbeddingEnergy(void)
{
	double rdata[17] = {
/*  1 */	0.0,
/*  2 */	0.0664348,
/*  3 */	0.0790068,
/*  4 */	0.0938295,
/*  5 */	0.111361,
/*  6 */	0.132086,
/*  7 */	0.156574,
/*  8 */	0.185503,
/*  9 */	0.21968,
/* 10 */	0.260071,
/* 11 */	0.307867,
/* 12 */	0.365323,
/* 13 */	0.435242,
/* 14 */	0.521086,
/* 15 */	0.627906,
/* 16 */	0.764184,
/* 17 */	0.9432,
	};
	
	double coeffdata[18*4] = {
/*  1 */	13.7590086032, 0.0, 0.0, 0.0,
/*  2 */	0.0, 13.7590086032, 0.0, 5.67698269175,
/*  3 */	0.0179989452968, 12.9462289287, 12.2342501826, -55.7077862349,
/*  4 */	-0.082894787807, 16.7773043944, -36.256172274, 148.875416841,
/*  5 */	0.0522393045964, 12.456675312, 9.79151102144, -14.7110344538,
/*  6 */	0.00600520789157, 13.7021893269, -1.39291253801, 18.7667984893,
/*  7 */	0.0248621935471, 13.2738999458, 1.84959569947, 10.5839668697,
/*  8 */	0.0303424824659, 13.1688959904, 2.5202313056, 9.15623615662,
/*  9 */	0.0192783958084, 13.3478272516, 1.55565714739, 10.8894963681,
/* 10 */	0.145696688406, 11.6214280623, 9.41436730519, -1.0350005328,
/* 11 */	-0.42320765032, 18.1839072248, -15.8190095716, 31.3066035492,
/* 12 */	0.521817948373, 8.97512916976, 14.0925587444, -1.07923214567,
/* 13 */	0.748417870515, 7.11431252151, 19.1861742021, -5.72681854372,
/* 14 */	0.621705672721, 7.98770307644, 17.1794979843, -4.18999196968,
/* 15 */	0.787828992214, 7.03129602755, 19.0149104678, -5.36408724229,
/* 16 */	-0.289324392858, 12.1777075834, 10.8187543558, -1.01303121786,
/* 17 */	6.31870309016, -13.7638052912, 44.7654498435, -15.8204191564,
/* 18 */	-6.95613165448, 28.4589592657, 0.0, 0.0};

	dArrayT		knots(17, rdata);
	dArray2DT	coeff(18, 4, coeffdata);

	fEmbeddingEnergy = new CubicSplineT(knots, coeff);
	if (!fEmbeddingEnergy) throw ExceptionT::kOutOfMemory;
}

/**********************************************************************
*
*  Glue functions
*
**********************************************************************/

/* private constructor */
VCPairPotential::VCPairPotential(void)
{
	/* parameters */
	fD		= 3.7760;
	falpha	= 1.4859;
	fR		= 2.1176;

	/* cut-off modifications */
	fq0		= 0.4204265806897907;
	fq1		=-0.06748375687440795;
}

/* I/O */
void VCPairPotential::Print(ostream& out) const
{
	out << " Parameters:\n\n";
	out << "     D = " << fD << '\n';
	out << " alpha = " << falpha << '\n';
	out << "     R = " << fR << '\n';

	out << " Cut-off modifications:\n\n";
	out << "    q1 = " << fq0 << '\n';
	out << "    q0 = " << fq1 << '\n';
}
	    	   	
void VCPairPotential::PrintName(ostream& out) const
{
	out << "    Voter-Chen pair potential\n";
}     	    	

/* returning values */
double VCPairPotential::Function(double x) const
{
	return( function(x) );	
}

double VCPairPotential::DFunction(double x) const
{
	return( Dfunction(x) );	
}

double VCPairPotential::DDFunction(double x) const
{
	return( DDfunction(x) );	
}

/* returning values in groups */
dArrayT& VCPairPotential::MapFunction(const dArrayT& in, dArrayT& out) const
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

dArrayT& VCPairPotential::MapDFunction(const dArrayT& in, dArrayT& out) const
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

dArrayT& VCPairPotential::MapDDFunction(const dArrayT& in, dArrayT& out) const
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
void VCPairPotential::SetAll(double x, dArrayT& data) const
{
	data[0] = function(x);
	data[1] = Dfunction(x);
	data[2] = DDfunction(x);
}   	

/**********************************************************************
* Private
**********************************************************************/

/* non-virtual function calls */
double VCPairPotential::function(double x) const
{
	if (x > kCutoffRadiusAl)
		return(0.0);
	else
		return( fD*(pow(1.0 - exp(falpha*(-x + fR)),2) - 1.0) + fq0 + fq1*x );
}

double VCPairPotential::Dfunction(double x) const
{
	if (x > kCutoffRadiusAl)
		return(0.0);
	else
	{
		double exp1 = exp(falpha*(-x + fR));
		return( 2.0*exp1*(1.0 - exp1)*falpha*fD + fq1 );		
	}
}

double VCPairPotential::DDfunction(double x) const
{
	if (x > kCutoffRadiusAl)
		return(0.0);
	else
	{
		double exp1 = exp(falpha*(x - fR));
		return( (2.0*(2.0 - exp1)*falpha*falpha*fD)/(exp1*exp1) );
	}
}	

/* private constructor */
VCElectronDensity::VCElectronDensity(void)
{
	/* parameters */
	fbeta = 3.3232;

	/* cut-off modifications */
	fq0 =-0.00380125678617228;
	fq1 = 0.0006334572870720792;
}

/* I/O */
void VCElectronDensity::Print(ostream& out) const
{
	out << " Parameters:\n\n";
	out << "  beta = " << fbeta << '\n';

	out << " Cut-off modifications:\n\n";
	out << "    q1 = " << fq0 << '\n';
	out << "    q0 = " << fq1 << '\n';
}
	    	   	
void VCElectronDensity::PrintName(ostream& out) const
{
	out << "    Voter-Chen electron density\n";
}     	    	

/* returning values */
double VCElectronDensity::Function(double x) const
{
	return( function(x) );	
}

double VCElectronDensity::DFunction(double x) const
{
	return( Dfunction(x) );	
}

double VCElectronDensity::DDFunction(double x) const
{
	return( DDfunction(x) );	
}

/* returning values in groups */
dArrayT& VCElectronDensity::MapFunction(const dArrayT& in, dArrayT& out) const
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

dArrayT& VCElectronDensity::MapDFunction(const dArrayT& in, dArrayT& out) const
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

dArrayT& VCElectronDensity::MapDDFunction(const dArrayT& in, dArrayT& out) const
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
void VCElectronDensity::SetAll(double x, dArrayT& data) const
{
	data[0] = function(x);
	data[1] = Dfunction(x);
	data[2] = DDfunction(x);
}   	

/**********************************************************************
* Private
**********************************************************************/

/* non-virtual function calls */
double VCElectronDensity::function(double x) const
{
	if (x > kCutoffRadiusAl)
		return(0.0);
	else
	{
		double exp1 = exp(-(x*fbeta));
	
		return( pow(x,6)*exp1*(512.0*exp1 + 1.0) + fq0 + x*fq1 );
	}
}

double VCElectronDensity::Dfunction(double x) const
{
	if (x > kCutoffRadiusAl)
		return(0.0);
	else
	{
		double exp1 = exp(-(x*fbeta));
				
		return( pow(x,5)*exp1*(6.0*(1.0 + 512.0*exp1) -
		                    fbeta*x*(1.0 + 1024.0*exp1)) + fq1
		       );
	}
}

double VCElectronDensity::DDfunction(double x) const
{
	if (x > kCutoffRadiusAl)
		return(0.0);
	else
	{
		double exp1 = exp(-(x*fbeta));
		
		return( pow(x,4)*exp1*(30.0*(1.0 +  512.0*exp1) -
					    12.0*x*fbeta*(1.0 + 1024.0*exp1) +
					 fbeta*fbeta*x*x*(1.0 + 2048.0*exp1))
			   );
	}
}	
