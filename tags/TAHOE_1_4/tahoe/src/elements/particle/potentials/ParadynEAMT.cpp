/* $Id: ParadynEAMT.cpp,v 1.9 2004-04-09 02:03:05 hspark Exp $ */
#include "ParadynEAMT.h"

#include "toolboxConstants.h"
#include "ifstreamT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"

using namespace Tahoe;

/* utility */
static inline int Min(int a, int b) { return (a < b) ? a : b; };
static inline double Min(double a, double b) { return (a < b) ? a : b; };

/* initialize static parameters */
int     ParadynEAMT::s_nr    = 0;
double  ParadynEAMT::s_f_inc = 1.0;
double* ParadynEAMT::s_Paircoeff = NULL;

int     ParadynEAMT::s_np    = 0;
double  ParadynEAMT::s_e_inc = 1.0;
double* ParadynEAMT::s_Embcoeff = NULL;

double* ParadynEAMT::s_ElecDenscoeff = NULL;

/* parameters */
const int knum_coeff = 9;

/* constructor */
ParadynEAMT::ParadynEAMT(const StringT& param_file):
  fParams(param_file),
  f_cut(0.0)
{
  const char caller[] = "ParadynEAMT::ParadynEAMT";
  /* try to open file */
  ifstreamT in(fParams);
  if (!in.is_open())
    ExceptionT::BadInputValue(caller, "error opening file: %s", fParams.Pointer());
  
  dArrayT tmp;	

  /* read comment line */
  fDescription.GetLineFromStream(in);
  
  /* lattice information */
  double mass;
  in >> fAtomicNumber >> mass 
     >> fLatticeParameter >> fStructure;

  /* Adjust mass like in interpolate.F of ParaDyn */
  double conmas = 1.0365e-4;
  mass *= conmas;
  
  /* read dimensions */
  int np, nr;
  double dp, dr;
  in >> np >> dp >> nr >> dr >> f_cut;
  if (np < 2   ||
      dp < 0.0 ||
      nr < 2   ||
      dr < 0.0 ||
      f_cut < 0.0) ExceptionT::BadInputValue(caller);
  
  /* Embedding Energy, frhoin in ParaDyn */
  tmp.Dimension(np);
  in >> tmp;
  /* compute spline coefficients for Embedded energy */
  ComputeCoefficients(tmp, dp, fEmbedCoeff);
  rho_inc = 1.0/dp;
  
  /* Pair Energy, zrin in ParaDyn 
     Note: It is only z at this point, not phi = z^2/r */
  tmp.Dimension(nr);
  in >> tmp;
  
  /* adjust units */
  tmp *= sqrt(27.2*0.529);
  
  f_inc = 1.0/dr;

  /* compute spline coefficients for z */
  ComputeCoefficients(tmp, dr, fPairCoeff);

  /* Electron Density, rhoin in ParaDyn, 
     assume that z and rho grids coincide */
  in >> tmp;
  /* compute spline coefficients for Electron Density  */
  ComputeCoefficients(tmp, dr, fElectronDensityCoeff);

  /* inherited */
  SetMass(mass);
  SetRange(f_cut);
}

/* write properties to output */

void ParadynEAMT::Write(ostream& out) const
{
  out << "Paradyn: " << fDescription << '\n';
  out << " Atomic number . . . . . . . . . . . . . . . . . = " << fAtomicNumber << '\n';
  out << " Lattice parameter . . . . . . . . . . . . . . . = " << fLatticeParameter << '\n';
  out << " Lattice structure . . . . . . . . . . . . . . . = " << fStructure << '\n';
  out << " Cut-off distance. . . . . . . . . . . . . . . . = " << f_cut << '\n';
  out << " # intervals in the electron density table . . . = " << fEmbedCoeff.MajorDim() << '\n';
  out << " # intervals in the potential table. . . . . . . = " << fPairCoeff.MajorDim() << '\n';
  out << " Interval size . . . . . . . . . . . . . . . . . = " << 1.0/f_inc << '\n';
}

/* return a pointer to the energy function */
ParadynEAMT::PairEnergyFunction ParadynEAMT::getPairEnergy(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_Paircoeff = fPairCoeff.Pointer();

  /* return function pointer */
  return ParadynEAMT::PairEnergy;
}

ParadynEAMT::EmbedEnergyFunction ParadynEAMT::getEmbedEnergy(void)
{
  /* copy my data to static */
  s_np       = fEmbedCoeff.MajorDim(); 
  s_e_inc    = rho_inc;
  s_Embcoeff = fEmbedCoeff.Pointer();

  /* return function pointer */
  return ParadynEAMT::EmbeddingEnergy;
}

ParadynEAMT::EDEnergyFunction ParadynEAMT::getElecDensEnergy(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_ElecDenscoeff = fElectronDensityCoeff.Pointer();
  
  /* return function pointer */
  return ParadynEAMT::ElecDensEnergy;
}

ParadynEAMT::PairForceFunction ParadynEAMT::getPairForce(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_Paircoeff = fPairCoeff.Pointer();

  /* return function pointer */
  return ParadynEAMT::PairForce;
}

ParadynEAMT::EmbedForceFunction ParadynEAMT::getEmbedForce(void)
{
  /* copy my data to static */
  s_np       = fEmbedCoeff.MajorDim(); 
  s_e_inc    = rho_inc;
  s_Embcoeff = fEmbedCoeff.Pointer();
  
  /* return function pointer */
  return ParadynEAMT::EmbeddingForce;
}

ParadynEAMT::EDForceFunction ParadynEAMT::getElecDensForce(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_ElecDenscoeff = fElectronDensityCoeff.Pointer();
  
  /* return function pointer */
  return ParadynEAMT::ElecDensForce;
}

ParadynEAMT::PairStiffnessFunction ParadynEAMT::getPairStiffness(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_Paircoeff = fPairCoeff.Pointer();

  /* return function pointer */
  return ParadynEAMT::PairStiffness;
}

ParadynEAMT::EmbedStiffnessFunction ParadynEAMT::getEmbedStiffness(void)
{
  /* copy my data to static */
  s_np       = fEmbedCoeff.MajorDim(); 
  s_e_inc    = rho_inc;
  s_Embcoeff = fEmbedCoeff.Pointer();
  
  /* return function pointer */
  return ParadynEAMT::EmbeddingStiffness;
}

ParadynEAMT::EDStiffnessFunction ParadynEAMT::getElecDensStiffness(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_ElecDenscoeff = fElectronDensityCoeff.Pointer();
  
  /* return function pointer */
  return ParadynEAMT::ElecDensStiffness;
}

/* return Paradyn-style coefficients table */
bool ParadynEAMT::getParadynTable(const double** coeff, double& dr, 
		int& row_size, int& num_rows) const
{
#pragma unused(coeff)
#pragma unused(dr)
#pragma unused(row_size)
#pragma unused(num_rows)
  return false;
}

/***********************************************************************
 * Private
 ***********************************************************************/

// z(r)
double ParadynEAMT::PairEnergy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  double z = EnergyAux(r_ab,s_nr,s_f_inc,s_Paircoeff);
  return z;
}

// F(rho)
double ParadynEAMT::EmbeddingEnergy(double rho_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)
  
  return EnergyAux(rho_ab,s_np,s_e_inc,s_Embcoeff);
}

// rho(r)
double ParadynEAMT::ElecDensEnergy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)
	
  return EnergyAux(r_ab,s_nr,s_f_inc,s_ElecDenscoeff);
}

// z'(r)
double ParadynEAMT::PairForce(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  double zp = ForceAux(r_ab,s_nr,s_f_inc,s_Paircoeff);
  return zp;
}

// F'(rho)
double ParadynEAMT::EmbeddingForce(double rho_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return ForceAux(rho_ab,s_np,s_e_inc,s_Embcoeff);
}

// rho(r)'
double ParadynEAMT::ElecDensForce(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return ForceAux(r_ab,s_nr,s_f_inc,s_ElecDenscoeff);
}

// z''(r)
double ParadynEAMT::PairStiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  double zpp = StiffnessAux(r_ab,s_nr,s_f_inc,s_Paircoeff);
  return zpp;
}

// F''(rho)
double ParadynEAMT::EmbeddingStiffness(double rho_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return StiffnessAux(rho_ab,s_np,s_e_inc,s_Embcoeff);
}

// rho''(r)
double ParadynEAMT::ElecDensStiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return StiffnessAux(r_ab,s_nr,s_f_inc,s_ElecDenscoeff);
}


/************************************************************************/
/* Compute energy, force and stiffness for a given atom i 
 * like in force.F */
double ParadynEAMT::EnergyAux(double r_ab,int n, double inc, double* coeff)
{
  double pp = r_ab*inc;
  int kk = int(pp);
  kk = Min(kk, n-2);
  pp -= kk;
  pp = Min(pp, 1.0);
  double* c = coeff + kk*knum_coeff;
  return c[0] + pp*(c[1] + pp*(c[2] + pp*c[3]));
}

double ParadynEAMT::ForceAux(double r_ab,int n, double inc, double* coeff)
{

  double pp = r_ab*inc;
  int kk = int(pp);
  kk = Min(kk, n-2);
  pp -= kk;
  pp = Min(pp, 1.0);
  double* c = coeff + kk*knum_coeff;
  return c[4] + pp*(c[5] + pp*c[6]);
}

double ParadynEAMT::StiffnessAux(double r_ab,int n, double inc, double* coeff)
{
  double pp = r_ab*inc;
  int kk = int(pp);
  kk = Min(kk, n-2);
  pp -= kk;
  pp = Min(pp, 1.0);
  double* c = coeff + kk*knum_coeff;
  return c[7] + pp*c[8];
}


/* compute the coefficients, like interpolation.F*/
void ParadynEAMT::ComputeCoefficients(const ArrayT<double>& f, double dx, dArray2DT& coeff)
{
  int nrar = f.Length();
	
  /* dimension */
  coeff.Dimension(nrar, knum_coeff);

  /* copy in function value */
  for (int j = 0; j < nrar; j++) coeff(j,0) = f[j];
  
  /* set function derivative at endpoints */
  coeff(0,1)      =      coeff(1,0)      - coeff(0,0);
  coeff(1,1)      = 0.5*(coeff(2,0)      - coeff(0,0));
  coeff(nrar-2,1) = 0.5*(coeff(nrar-1,0) - coeff(nrar-3,0));
  // Syl: coeff(nrar-1,1) = 0.0;
  coeff(nrar-1,1) = coeff(nrar-1,0) - coeff(nrar-2,0);
  
  /* derivative approximation through the middle */
  for (int j = 2; j < nrar-2; j++)
    coeff(j,1) = ((coeff(j-2,0) - coeff(j+2,0)) + 
              8.0*(coeff(j+1,0) - coeff(j-1,0)))/12.0;
  
  /* higher order coefficients */
  for (int j = 0; j < nrar-1; j++)
    {
      coeff(j,2) = 3.0*(coeff(j+1,0) - coeff(j,0)) - 
                   2.0*coeff(j,1) - coeff(j+1,1);
      coeff(j,3) = coeff(j,1) + coeff(j+1,1) - 
                   2.0*(coeff(j+1,0) - coeff(j,0));
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
