/* $Id: ParadynT.cpp,v 1.1 2003-04-03 01:13:00 saubry Exp $ */
#include "ParadynT.h"

#include "toolboxConstants.h"
#include "ifstreamT.h"
#include "dArrayT.h"
#include "AutoArrayT.h"

using namespace Tahoe;

/* utility */
static inline int Min(int a, int b) { return (a < b) ? a : b; };
static inline double Min(double a, double b) { return (a < b) ? a : b; };

/* initialize static parameters */
int     ParadynT::s_nr    = 0;
double  ParadynT::s_f_inc = 1.0;
double* ParadynT::s_Paircoeff = NULL;

int     ParadynT::s_np    = 0;
double  ParadynT::s_e_inc = 1.0;
double* ParadynT::s_Embcoeff = NULL;

double* ParadynT::s_ElecDenscoeff = NULL;

/* parameters */
const int knum_coeff = 9;

/* constructor */
ParadynT::ParadynT(const StringT& param_file):
  fParams(param_file),
  f_cut(0.0)
{
  const char caller[] = "ParadynT::ParadynT";
  
  /* try to open file */
  ifstreamT in(fParams);
  if (!in.is_open())
    ExceptionT::BadInputValue(caller, "error opening file: %s", fParams.Pointer());
  
  /* read comment line */
  fDescription.GetLineFromStream(in);
  
  /* lattice information */
  double mass;
  in >> fAtomicNumber >> mass 
     >> fLatticeParameter >> fStructure;
  
  /* read dimensions */
  int np, nr;
  double dp, dr;
  in >> np >> dp >> nr >> dr >> f_cut;
  if (np < 2   ||
      dp < 0.0 ||
      nr < 2   ||
      dr < 0.0 ||
      f_cut < 0.0) ExceptionT::BadInputValue(caller);
  
  /* read tables of F(rho),z(r) and rho(r) */
  dArrayT tmp;	

  /* Embedding Energy, frhoin in ParaDyn */
  tmp.Dimension(np);
  in >> tmp;
  /* compute spline coefficients for Embedded energy */
  ComputeCoefficients(tmp, dp, fEmbedCoeff);
  rho_inc = 1.0/dp;
    
  /* Pair Energy, zrin in ParaDyn  */
  tmp.Dimension(nr);
  in >> tmp;
  tmp *= sqrt(27.2*0.529);

  /* convert z to phi */
  dArrayT phi;
  phi.Dimension(nr);
  f_inc = 1.0/dr;

  phi[0] = tmp[0]*tmp[0]*f_inc;
  double r = dr;
  for (int j = 1; j < nr; j++)
    {
      phi[j] = tmp[j]*tmp[j]/r;      
      r += dr;
    }
  /* compute spline coefficients for Pair energy: phi */
  ComputeCoefficients(phi, dr, fPairCoeff);


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

void ParadynT::Write(ostream& out) const
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
ParadynT::PairEnergyFunction ParadynT::getPairEnergy(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_Paircoeff = fPairCoeff.Pointer();

  /* return function pointer */
  return ParadynT::PairEnergy;
}

ParadynT::EmbedEnergyFunction ParadynT::getEmbedEnergy(void)
{
  /* copy my data to static */
  s_np       = fEmbedCoeff.MajorDim(); 
  s_e_inc    = rho_inc;
  s_Embcoeff = fEmbedCoeff.Pointer();

  /* return function pointer */
  return ParadynT::EmbeddingEnergy;
}

ParadynT::EDEnergyFunction ParadynT::getElecDensEnergy(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_ElecDenscoeff = fElectronDensityCoeff.Pointer();
  
  /* return function pointer */
  return ParadynT::ElecDensEnergy;
}

ParadynT::PairForceFunction ParadynT::getPairForce(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_Paircoeff = fPairCoeff.Pointer();

  /* return function pointer */
  return ParadynT::PairForce;
}

ParadynT::EmbedForceFunction ParadynT::getEmbedForce(void)
{
  /* copy my data to static */
  s_np       = fEmbedCoeff.MajorDim(); 
  s_e_inc    = rho_inc;
  s_Embcoeff = fEmbedCoeff.Pointer();
  
  /* return function pointer */
  return ParadynT::EmbeddingForce;
}

ParadynT::EDForceFunction ParadynT::getElecDensForce(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_ElecDenscoeff = fElectronDensityCoeff.Pointer();
  
  /* return function pointer */
  return ParadynT::ElecDensForce;
}

ParadynT::PairStiffnessFunction ParadynT::getPairStiffness(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_Paircoeff = fPairCoeff.Pointer();

  /* return function pointer */
  return ParadynT::PairStiffness;
}

ParadynT::EmbedStiffnessFunction ParadynT::getEmbedStiffness(void)
{
  /* copy my data to static */
  s_np       = fEmbedCoeff.MajorDim(); 
  s_e_inc    = rho_inc;
  s_Embcoeff = fEmbedCoeff.Pointer();
  
  /* return function pointer */
  return ParadynT::EmbeddingStiffness;
}

ParadynT::EDStiffnessFunction ParadynT::getElecDensStiffness(void)
{
  /* copy my data to static */
  s_nr        = fPairCoeff.MajorDim(); 
  s_f_inc     = f_inc;
  s_ElecDenscoeff = fElectronDensityCoeff.Pointer();
  
  /* return function pointer */
  return ParadynT::ElecDensStiffness;
}

/* return Paradyn-style coefficients table */
bool ParadynT::getParadynTable(const double** coeff, double& dr, int& row_size, int& num_rows) const
{
	*coeff = fPairCoeff.Pointer();
	dr = f_inc;
	row_size = 9;
	num_rows = fPairCoeff.MajorDim();
	return true;
}

/***********************************************************************
 * Private
 ***********************************************************************/

// phi
double ParadynT::PairEnergy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return  EnergyAux(r_ab,s_nr,s_f_inc,s_Paircoeff);
}

// F(rho)
double ParadynT::EmbeddingEnergy(double rho_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return EnergyAux(rho_ab,s_np,s_e_inc,s_Embcoeff);
}

// rho
double ParadynT::ElecDensEnergy(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return EnergyAux(r_ab,s_nr,s_f_inc,s_ElecDenscoeff);
}

// f'
double ParadynT::PairForce(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return ForceAux(r_ab,s_nr,s_f_inc,s_Paircoeff);
}

// F'(rho)
double ParadynT::EmbeddingForce(double rho_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return ForceAux(rho_ab,s_np,s_e_inc,s_Embcoeff);
}

// rho'
double ParadynT::ElecDensForce(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return ForceAux(r_ab,s_nr,s_f_inc,s_ElecDenscoeff);
}

// f''
double ParadynT::PairStiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return StiffnessAux(r_ab,s_nr,s_f_inc,s_Paircoeff);
}

// F''(rho)
double ParadynT::EmbeddingStiffness(double rho_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return StiffnessAux(rho_ab,s_np,s_e_inc,s_Embcoeff);
}

// rho''
double ParadynT::ElecDensStiffness(double r_ab, double* data_a, double* data_b)
{
#pragma unused(data_a)
#pragma unused(data_b)

  return StiffnessAux(r_ab,s_nr,s_f_inc,s_ElecDenscoeff);
}


/************************************************************************/
/* Compute energy, force and stiffness for a given atom i 
 * like in force.F */
double ParadynT::EnergyAux(double r_ab,int n, double inc, double* coeff)
{
  double pp = r_ab*inc;
  int kk = int(pp);
  kk = Min(kk, n-2);
  pp -= kk;
  pp = Min(pp, 1.0);
  double* c = coeff + kk*knum_coeff;
  return c[0] + pp*(c[1] + pp*(c[2] + pp*c[3]));
}

double ParadynT::ForceAux(double r_ab,int n, double inc, double* coeff)
{

  double pp = r_ab*inc;
  int kk = int(pp);
  kk = Min(kk, n-2);
  pp -= kk;
  pp = Min(pp, 1.0);
  double* c = coeff + kk*knum_coeff;
  return c[4] + pp*(c[5] + pp*c[6]);
}

double ParadynT::StiffnessAux(double r_ab,int n, double inc, double* coeff)
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
void ParadynT::ComputeCoefficients(const ArrayT<double>& f, double dx, dArray2DT& coeff)
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
  coeff(nrar-1,1) = 0.0;
  
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
