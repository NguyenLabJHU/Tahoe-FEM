/* $Id: ParadynPairT.cpp,v 1.6 2003-01-27 07:00:26 paklein Exp $ */
#include "ParadynPairT.h"
#include "toolboxConstants.h"
#include "ifstreamT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"

using namespace Tahoe;

/* utility */
static inline int Min(int a, int b) { return (a < b) ? a : b; };
static inline double Min(double a, double b) { return (a < b) ? a : b; };

/* initialize static parameters */
int     ParadynPairT::s_nr    = 0;
double  ParadynPairT::s_1bydr = 1.0;
double* ParadynPairT::s_coeff = NULL;

/* parameters */
const int knum_coeff = 9;

/* constructor */
ParadynPairT::ParadynPairT(const StringT& param_file):
	fParams(param_file),
	f_cut(0.0)
{
	const char caller[] = "ParadynPairT::ParadynPairT";

	/* try to open file */
	ifstreamT in(fParams);
	if (!in.is_open())
		ExceptionT::BadInputValue(caller, "error opening file: %s", fParams.Pointer());

	/* read comment line */
	fDescription.GetLineFromStream(in);

	/* lattice information */
	double mass;
	in >> fAtomicNumber >> mass >> fLatticeParameter >> fStructure;
	
	/* table dimensions */
	int np, nr;
	double dp, dr;
	in >> np >> dp >> nr >> dr >> f_cut;
	if (np < 2   ||
	    dp < 0.0 ||
	    nr < 2   ||
	    dr < 0.0 ||
	 f_cut < 0.0) ExceptionT::BadInputValue(caller);
	
	/* embedding energy - not used */
	dArrayT tmp(np);
	in >> tmp;
	
	/* phi function */
	tmp.Dimension(nr);
	in >> tmp;

	/* compute spline coefficients */
	ComputeCoefficients(tmp, dr, fCoefficients);
	f_1bydr = 1.0/dr;

	/* inherited */
	SetMass(mass);
	SetRange(f_cut);
}

/* write properties to output */
void ParadynPairT::Write(ostream& out) const
{
	/* inherited */
	PairPropertyT::Write(out);

	out << "Paradyn: " << fDescription << '\n';
	out << " Atomic number . . . . . . . . . . . . . . . . . = " << fAtomicNumber << '\n';
	out << " Lattice parameter . . . . . . . . . . . . . . . = " << fLatticeParameter << '\n';
	out << " Lattice structure . . . . . . . . . . . . . . . = " << fStructure << '\n';
	out << " Cut-off distance. . . . . . . . . . . . . . . . = " << f_cut << '\n';
	out << " Number of intervals in the potential table. . . = " << fCoefficients.MajorDim() << '\n';
	out << " Interval size . . . . . . . . . . . . . . . . . = " << 1.0/f_1bydr << '\n';
}

/* return a pointer to the energy function */
PairPropertyT::EnergyFunction ParadynPairT::getEnergyFunction(void)
{
	/* copy my data to static */
	s_nr    = fCoefficients.MajorDim(); 
	s_1bydr = f_1bydr;
	s_coeff = fCoefficients.Pointer();

	/* return function pointer */
	return ParadynPairT::Energy;
}

PairPropertyT::ForceFunction ParadynPairT::getForceFunction(void)
{
	/* copy my data to static */
	s_nr    = fCoefficients.MajorDim(); 
	s_1bydr = f_1bydr;
	s_coeff = fCoefficients.Pointer();

	/* return function pointer */
	return ParadynPairT::Force;
}

PairPropertyT::StiffnessFunction ParadynPairT::getStiffnessFunction(void)
{
	/* copy my data to static */
	s_nr    = fCoefficients.MajorDim(); 
	s_1bydr = f_1bydr;
	s_coeff = fCoefficients.Pointer();

	/* return function pointer */
	return ParadynPairT::Stiffness;
}

/* return Paradyn-style coefficients table */
bool ParadynPairT::getParadynTable(const double** coeff, double& dr, int& row_size, int& num_rows) const
{
	*coeff = fCoefficients.Pointer();
	dr = f_1bydr;
	row_size = 9;
	num_rows = fCoefficients.MajorDim();
	return true;
}

/***********************************************************************
 * Private
 ***********************************************************************/

double ParadynPairT::Energy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double pp = r_ab*s_1bydr;
	int kk = int(pp);
	kk = Min(kk, s_nr-2);
	pp -= kk;
	pp = Min(pp, 1.0);
	double* c = s_coeff + kk*knum_coeff;
	return c[0] + pp*(c[1] + pp*(c[2] + pp*c[3]));
}

double ParadynPairT::Force(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double pp = r_ab*s_1bydr;
	int kk = int(pp);
	kk = Min(kk, s_nr-2);
	pp -= kk;
	pp = Min(pp, 1.0);
	double* c = s_coeff + kk*knum_coeff;
	return c[4] + pp*(c[5] + pp*c[6]);
}

double ParadynPairT::Stiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

	double pp = r_ab*s_1bydr;
	int kk = int(pp);
	kk = Min(kk, s_nr-2);
	pp -= kk;
	pp = Min(pp, 1.0);
	double* c = s_coeff + kk*knum_coeff;
	return c[7] + pp*c[8];
}

/* compute the coefficients */
void ParadynPairT::ComputeCoefficients(const ArrayT<double>& f, double dx, dArray2DT& coeff)
{
	int nrar = f.Length();
	
	/* dimension */
	coeff.Dimension(nrar, knum_coeff);

	/* copy in function value */
	for (int j = 0; j < nrar; j++)
		coeff(j,0) = f[j];

	/* set function derivative at endpoints */
	coeff(0,1) = coeff(1,0) - coeff(0,0);
	coeff(1,1) = 0.5*(coeff(2,0) - coeff(0,0));
	coeff(nrar-2,1) = 0.5*(coeff(nrar-1,0) - coeff(nrar-3,0));
	coeff(nrar-1,1) = 0.0;

	/* derivative approximation through the middle */
	for (int j = 2; j < nrar-2; j++)
		coeff(j,1) = ((coeff(j-2,0) - coeff(j+2,0)) + 8.0*(coeff(j+1,0) - coeff(j-1,0)))/12.0;

	/* higher order coefficients */
	for (int j = 0; j < nrar-1; j++)
	{
		coeff(j,2) = 3.0*(coeff(j+1,0) - coeff(j,0)) - 2.0*coeff(j,1) - coeff(j+1,1);
		coeff(j,3) = coeff(j,1) + coeff(j+1,1) - 2.0*(coeff(j+1,0) - coeff(j,0));
	}
	coeff(nrar-1,2) = 0.0;
	coeff(nrar-1,3) = 0.0;
    
    /* coefficients for derivatives */
	for (int j = 0; j < nrar; j++)
	{
		/* for first derivative */
		coeff(j,4) = coeff(j,1)/dx;
		coeff(j,5) = 2.0*coeff(j,2)/dx;
		coeff(j,6) = 3.0*coeff(j,3)/dx;

		/* for second derivatives */
		coeff(j,7) = coeff(j,5)/dx;
		coeff(j,8) = 2.0*coeff(j,6)/dx;
	}
}
