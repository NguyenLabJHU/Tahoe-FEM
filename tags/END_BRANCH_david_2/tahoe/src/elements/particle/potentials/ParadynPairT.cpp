/* $Id: ParadynPairT.cpp,v 1.9 2004-07-15 08:29:49 paklein Exp $ */
#include "ParadynPairT.h"

#include "ifstreamT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"
#include "BasicSupportT.h"

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
ParadynPairT::ParadynPairT(const BasicSupportT* support, const StringT& param_file):
	fSupport(support),
	f_cut(0.0)
{
	SetName("Paradyn_pair");

	/* initialize potential from stream */
	ReadParameters(param_file);
}

ParadynPairT::ParadynPairT(const BasicSupportT* support):
	fSupport(support)
{
	SetName("Paradyn_pair");
}

#if 0
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
#endif

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

/* describe the parameters needed by the interface */
void ParadynPairT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	PairPropertyT::DefineParameters(list);

	/* give "mass" default value */
	ParameterT& mass = list.GetParameter("mass");
	mass.SetDefault(1.0);

	/* parameter file path (relative to input file) */
	list.AddParameter(ParameterT::Word, "parameter_file");
}

/* accept parameter list */
void ParadynPairT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	PairPropertyT::TakeParameterList(list);

	/* convert file path to standard form */
	StringT file = list.GetParameter("parameter_file");
	file.ToNativePathName();
	
	/* prepend path from input file */
	if (!fSupport) ExceptionT::GeneralFail("ParadynPairT::TakeParameterList",
		"pointer to BasicSupportT not set");
	StringT path;
	path.FilePath(fSupport->InputFile());
	file.Prepend(path);

	/* read parameters */
	ReadParameters(file);
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
  // Syl: coeff(nrar-1,1) = 0.0;
  coeff(nrar-1,1) = coeff(nrar-1,0) - coeff(nrar-2,0);
  
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

/* initialize potential from stream */
void ParadynPairT::ReadParameters(const StringT& param_file)
{
	const char caller[] = "ParadynPairT::ReadParameters";

	/* try to open file */
	ifstreamT in(param_file);
	if (!in.is_open())
		ExceptionT::BadInputValue(caller, "error opening file: %s", param_file.Pointer());

	/* read comment line */
	fDescription.GetLineFromStream(in);

	/* lattice information */
	double mass;
	in >> fAtomicNumber >> mass >> fLatticeParameter >> fStructure;
	
	/* Adjust mass like in interpolate_pair.F of ParaDyn */
	double conmas = 1.0365e-4;
	mass *= conmas;

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

	/* set parmeters */
	SetMass(mass);
	SetRange(f_cut);
}
