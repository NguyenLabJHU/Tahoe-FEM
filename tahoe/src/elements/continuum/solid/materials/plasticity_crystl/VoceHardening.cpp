/*
  File: VoceHardening.cpp
*/

#include "VoceHardening.h"
#include "PolyCrystalMatT.h"
#include "NLCSolver.h"
#include "Utils.h"


// number of material properties and initial hard values

using namespace Tahoe;

const int kNumMatProp    = 5;
const int kNumInitValues = 1;
const int kNumInternal   = 2;

// codes for computing hardening law qnts
const int kFunc  = 0;  // HardFunc 
const int kdFunc = 1;  // d(HardFunc)/d(Crss)

// codes to access fInternal
const int kShearRate_n = 0;
const int kShearRate   = 1;

// some limits for total shear rate
const double SHR_MIN = 0.;
const double SHR_MAX = 1.0e6;

VoceHardening::VoceHardening(PolyCrystalMatT& poly) :
  SlipHardening(poly),
  fInternal (kNumInternal),
  farray    (poly.NumSlip())
{
  // input file
  ifstreamT& in = poly.Input_x();

  // number of hardening variables
  // fNumHardVar = poly.NumSlip();
  fNumHardVar = 1;

  // allocate space for
  // ... material constants and initial hard values
  fMatProp.Dimension(kNumMatProp);
  fInitHardValues.Dimension(kNumInitValues);

  // ... hardening variables
  fTauIso.Dimension(fNumHardVar);
  fTauIso_n.Dimension(fNumHardVar);
  ftauiso_save.Dimension(fNumHardVar);

  // input material properties for hardening law
  in >> fMatProp[0];   // h0
  in >> fMatProp[1];   // g0
  in >> fMatProp[2];   // gs0
  in >> fMatProp[3];   // gams0
  in >> fMatProp[4];   // xms

  // input initial slip system hardness (same for all) 
  in >> fInitHardValues[0];

  // set hardening solver (in base class)
  SetHardeningSolver(in, fNumHardVar);
}

VoceHardening::~VoceHardening() {}

void VoceHardening::InitializeHardVariables() 
{ 
  // hardening stress
  fTauIso_n = fInitHardValues[0]; 
  fTauIso   = fInitHardValues[0]; 

  // internal variables
  fInternal = 0.0;
}

double VoceHardening::Magnitude() const { return fTauIso.Magnitude(); }
const int VoceHardening::NumberOfVariables() const 
{ 
  return 2*fNumHardVar + kNumInternal; 
}

void VoceHardening::UpdateHistory() 
{ 
  fTauIso_n = fTauIso; 
  fInternal[kShearRate_n] = fInternal[kShearRate];
}

void VoceHardening::ResetHistory() { 
  fTauIso = fTauIso_n; 
  fInternal[kShearRate] = fInternal[kShearRate_n];
}

void VoceHardening::LoadHardData(int dim, int dex, dArrayT& d_array)
{
  // recover hardening variables for current element/IP/grain
  fTauIso_n.Set (fNumHardVar,  &d_array[dex += dim        ]);
  fTauIso.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
  fInternal.Set (kNumInternal, &d_array[dex += fNumHardVar]);
}

void VoceHardening::ExplicitUpdateHard()
{
  // preliminary computations
  InternalHardQnts();

  // forward Euler estimate for crss
  double c   = fdt * fMatProp[0];
  double g_s = fTauIsoSat - fMatProp[1];
  double g_n, g;
  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i] - fMatProp[1];

     if ( (g_n/g_s) <= 1.0 )
        g = g_n + c * (1.0 - g_n/g_s) * fInternal[kShearRate];
     else
        g = g_n;

     fTauIso[i] = g + fMatProp[1];
  }

  // norm of explicit estimate
  fNormHard0 = fTauIso.Magnitude();
}

void VoceHardening::ImplicitUpdateHard()   // called from Algorithm 1
{
  // preliminary computations
  InternalHardQnts();

  // generalized mid-point approximation for crss (theta=0.5)
  double c   = 0.5 * fdt * fMatProp[0];
  double g_s = fTauIsoSat - fMatProp[1];
  double g_n, g;

  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i] - fMatProp[1];
     if ( (g_n/g_s) <= 1.0 )
     {
        g = g_n + c*( (1.0-g_n/g_s)*fInternal[kShearRate_n] + fInternal[kShearRate] );
        g /= ( 1.0 + c*fInternal[kShearRate]/g_s );
     }
     else
     {
        g = g_n;
     }

     fTauIso[i] = g + fMatProp[1];
  }

  // norm of implicit estimate
  fNormHard = fTauIso.Magnitude();
}

void VoceHardening::ImplicitSolveHard()   // called from Algorithm 2
{
  int ierr = 0;
  
  // preliminary computations
  InternalHardQnts();

  // backward Euler method: call solver method of NLCSolver class
  //fSolver->Solve(fSolverPtr, fTauIso, ierr);

  //if (ierr == 1)
  //  throwRunTimeError("VoceHardening::SolveImplicitHard: Convergence problems");

  // backward Euler approximation for crss
  double c   = fdt * fMatProp[0];
  double g_s = fTauIsoSat - fMatProp[1];
  double g_n, g;

  for (int i = 0; i < fNumHardVar; i++)
  {
     g_n = fTauIso_n[i] - fMatProp[1];

     if ( (g_n/g_s) <= 1.0 )
     {
        g = g_n + c*fInternal[kShearRate];
        g /= ( 1.0 + c*fInternal[kShearRate]/g_s );
     }
     else
     {
        g = g_n;
     }

     fTauIso[i] = g + fMatProp[1];
  }

  // norm of implicit solution
  fNormHard = fTauIso.Magnitude();
}

void VoceHardening::FormRHS(const dArrayT& tauIso, dArrayT& rhs)
{
  // form residual
  for (int i = 0; i < fNumHardVar; i++)
    rhs[i] = tauIso[i] - fTauIso_n[i] - fdt * HardeningLaw(tauIso[i], kFunc);
}

void VoceHardening::FormLHS(const dArrayT& tauIso, dMatrixT& lhs)
{
  // form jacobian
  lhs = 0.0;
  for (int i = 0; i < fNumHardVar; i++)
    lhs(i,i) = 1. - fdt * HardeningLaw(tauIso[i], kdFunc);
}

bool VoceHardening::Converged(double toler)
{
  // check convergence on hardening variables
  bool test = ( fabs(fNormHard-fNormHard0) < toler*fInitHardValues[0] );

  // if did not converge, reset norm
  if (!test) fNormHard0 = fNormHard;

  return test;
} 

void VoceHardening::SaveCurrentSolution() { ftauiso_save = fTauIso; }
void VoceHardening::RestoreSavedSolution() { fTauIso = ftauiso_save; }

double VoceHardening::HardeningModulus() const
{
  // for now return h0
  return  fMatProp[0];
}

const double VoceHardening::IsoHardeningStress(int is) const
{
  // one isotropic hardening variable per grain
  #pragma unused(is)
  return fTauIso[0];
}

const dArrayT& VoceHardening::ComputeHardQnts()
{
  // preliminary computations
  InternalHardQnts();

  // compute dHardLaw/dTauIso
  double dHdIso = fdt * HardeningLaw(fTauIso[0], kdFunc);

  // compute factor of dHardLaw/dDGamma
  double dHdDGam = HardeningLaw(fTauIso[0], kFunc) / (fInternal[kShearRate]*fdt);

  // compute array (1-dHardLaw/dTauIso)^(-1)*dHardLaw/dDGamma
  for (int i = 0; i < fDGamma.Length(); i++)
    {
      farray[i] = fdt * dHdDGam * fabs(fDGamma[i]) / fDGamma[i];
      farray[i] *= 1./ (1. - dHdIso);
    }
  return farray;
}

#if 0
  // print hardening parameters
  out << "       Hardening rate (h0) . . . . . . . . . . . = " << fMatProp[0] << "\n";
  out << "       Initial hardness (g0) . . . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       Saturation hardness (gs0) . . . . . . . . = " << fMatProp[2] << "\n";
  out << "       Saturation shear rate (gams0) . . . . . . = " << fMatProp[3] << "\n";
  out << "       Saturation exponent (xms) . . . . . . . . = " << fMatProp[4] << "\n";
#endif

/* PRIVATE MEMBERS FUNCTIONS */

void VoceHardening::InternalHardQnts()
{
  // accumulated shear rate
  fInternal[kShearRate] = 0.;
  for (int i = 0; i < fDGamma.Length(); i++) fInternal[kShearRate] += fabs(fDGamma[i]);
  fInternal[kShearRate] /= fdt; 

  // check limits on fShearRate
  if (fInternal[kShearRate] <= SHR_MIN) fInternal[kShearRate] = SHR_MIN;
  if (fInternal[kShearRate] >= SHR_MAX) fInternal[kShearRate] = SHR_MAX;

  // hardening saturation level
  fTauIsoSat = fMatProp[2] * pow(fInternal[kShearRate]/fMatProp[3], fMatProp[4]);
}
	
const double VoceHardening::HardeningLaw(double tauIso, int kcode)
{
  double hard = fMatProp[0] / (fTauIsoSat - fMatProp[1]) * fInternal[kShearRate];

  switch(kcode)
    {
    case kFunc:
      hard *= (fTauIsoSat - tauIso);
      break;
    case kdFunc:
      hard *= -1.0;
      break;
    default:
      throwRunTimeError("VoceHardening::HardeningLaw: Bad kcode");
    }

  return hard;
}
