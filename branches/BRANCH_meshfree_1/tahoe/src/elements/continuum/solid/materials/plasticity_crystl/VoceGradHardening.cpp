/*
  File: VoceGradHardening.cpp
*/

#include "VoceGradHardening.h"
#include "PolyCrystalMatT.h"
#include "NLCSolver.h"
#include "Utils.h"
#include "ifstreamT.h"

// number of material properties for hard model 
const int kNumMatProp = 6;
const int kNumInitValues = 1;

// codes for computing hardening law qnts
const int kFunc  = 0;  // HardFunc 
const int kdFunc = 1;  // d(HardFunc)/d(tauS)

VoceGradHardening::VoceGradHardening(PolyCrystalMatT& poly):
  SlipHardening(poly)
{
  // input file
  ifstreamT& in = poly.Input_x();

  // number of hardening variables
  fNumHardVar = poly.NumSlip();

  // allocate space for
  // ... material constants and initial hard values
  fMatProp.Allocate(kNumMatProp);
  fInitHardValues.Allocate(kNumInitValues);

  // ... hardening variables
  fTauIso.Allocate(fNumHardVar);
  fTauKin.Allocate(fNumHardVar);
  fTauIso_n.Allocate(fNumHardVar);

  // input material properties for hardening law
  in >> fMatProp[0];   // h0 
  in >> fMatProp[1];   // tauS0
  in >> fMatProp[2];   // tauSsat
  in >> fMatProp[3];   // c_s
  in >> fMatProp[4];   // c_x
  in >> fMatProp[5];   // burger's vector b

  // input initial value of statistically stored dislocation
  in >> fInitHardValues[0];

  // set hardening solver (in base class)
  SetHardeningSolver(in, fNumHardVar);
}

VoceGradHardening::~VoceGradHardening() {}

void VoceGradHardening::InitializeHardVariables() 
{ 
  // initialize hardening variables at t_0
  fTauIso_n = fInitHardValues[0];
  ResetHistory();
}

double VoceGradHardening::Magnitude() const { return fTauIso.Magnitude(); }
const int VoceGradHardening::NumberOfVariables() const { return 2*fNumHardVar; }

void VoceGradHardening::UpdateHistory()
{
  // update hardening variables
  fTauIso_n = fTauIso;
  //fTauKin_n = fTauKin;  // keep track of fXe (in GradCrystalPlast)
}

void VoceGradHardening::ResetHistory()
{
  // reset hardening variables
  fTauIso = fTauIso_n;
  //fTauKin = fTauKin_n;  // keep track of fXe (in GradCrystalPlast)
}

void VoceGradHardening::LoadHardData(int dim, int dex, dArrayT& d_array)
{
  // recover hardening variables for current element/IP/grain
  fTauIso_n.Set (fNumHardVar,  &d_array[dex += dim        ]);
  fTauIso.Set   (fNumHardVar,  &d_array[dex += fNumHardVar]);
}

void VoceGradHardening::ExplicitUpdateHard()
{
  // preliminary computations
  InternalHardQnts();

  // forward Euler estimate for fTauS
  for (int i = 0; i < fNumHardVar; i++)
    fTauIso[i] = fTauIso_n[i] + fdt * HardeningLaw(fTauIso_n[i], kFunc);
}

void VoceGradHardening::ImplicitUpdateHard()
{
  // preliminary computations
  InternalHardQnts();

  // forward Euler estimate for fTauS
  for (int i = 0; i < fNumHardVar; i++)
    fTauIso[i] = fTauIso_n[i] + fdt * HardeningLaw(fTauIso[i], kFunc);
}

void VoceGradHardening::ImplicitSolveHard()
{
  int ierr = 0;
  
  // preliminary computations
  InternalHardQnts();

  // compute DDss
  fSolver->Solve(fSolverPtr, fTauIso, ierr);

  if (ierr == 1)
    throwRunTimeError("VoceGradHardening::SolveImplicitHard: Convergence problems");
}

void VoceGradHardening::FormRHS(const dArrayT& tauIso, dArrayT& rhs)
{
  // form residual
  for (int i = 0; i < fNumHardVar; i++)
    rhs[i] = tauIso[i] - fTauIso_n[i] - fdt * HardeningLaw(tauIso[i], kFunc);
}

void VoceGradHardening::FormLHS(const dArrayT& tauIso, dMatrixT& lhs)
{
  // form Jacobian
  lhs = 0.0;
  for (int i = 0; i < fNumHardVar; i++)
    lhs(i,i) = 1. - fdt * HardeningLaw(tauIso[i], kdFunc);
}

bool VoceGradHardening::Converged(double toler)
{
  #pragma unused(toler)
  return false;
}

double VoceGradHardening::HardeningModulus() const
{
  // temporary value
  return  fMatProp[0];
}

void VoceGradHardening::Print(ostream& out) const
{
  // print hardening parameters
  out << "       Hardening rate (h0) . . . . . . . . . . . = " << fMatProp[0] << "\n";
  out << "       Initial hardness (tauS0). . . . . . . . . = " << fMatProp[1] << "\n";
  out << "       Saturation hardness (tauSsat) . . . . . . = " << fMatProp[2] << "\n";
  out << "       c_s . . . . . . . . . . . . . . . . . . . = " << fMatProp[3] << "\n";
  out << "       c_x . . . . . . . . . . . . . . . . . . . = " << fMatProp[4] << "\n";
  out << "       burger's vector (b) . . . . . . . . . . . = " << fMatProp[5] << "\n";

  // print hardening solver control data
  fSolver->Print(out);
}

void VoceGradHardening::PrintName(ostream& out) const
{
  // print model name
  out << "    VoceGrad's slip hardening law\n";
}

/* PRIVATE MEMBER FUNCTIONS */

void VoceGradHardening::InternalHardQnts()
{
  // accumulated shear rate and work rate
  fShearRate = 0.;
  fWorkRate  = 0.;
  for (int i = 0; i < fDGamma.Length(); i++)
    {
      fShearRate += fabs(fDGamma[i]);
      fWorkRate  += fabs(fTauKin[i])*fabs(fDGamma[i]);
    }
  fShearRate /= fdt;
  fWorkRate  /= fdt;

  // saturation value of hardness
  fTauIsoSat = fMatProp[2];
}

const double VoceGradHardening::HardeningLaw(double tauIso, int kcode)
{
  double coeff = 50.e3*fMatProp[3]*fMatProp[3]/(2.*fMatProp[4]); // mu*c_s^2/(2*c_x)
  double diff  = tauIso-fMatProp[1];
  double hard;

  switch(kcode)
    {
    case kFunc:
      hard = fMatProp[0]*(fTauIsoSat-tauIso)/(fTauIsoSat-fMatProp[1])*fShearRate;
      if (diff > 0.) hard += coeff/diff*fWorkRate;
      break;

    case kdFunc:
      hard = -fMatProp[0]/(fTauIsoSat-fMatProp[1])*fShearRate;
      if (diff > 0.) hard -= coeff/(diff*diff)*fWorkRate;
      break;

    default:
      throwRunTimeError("VoceGradHardening::HardeningLaw: Bad kcode");
    }

  return hard;
}
