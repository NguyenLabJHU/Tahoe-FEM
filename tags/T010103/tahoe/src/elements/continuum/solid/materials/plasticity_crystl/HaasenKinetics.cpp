/*
  File: HaasenKinetics.cpp
*/

#include "HaasenKinetics.h"
#include "PolyCrystalMatT.h"
#include "Utils.h"


using namespace Tahoe;

const int kNumMatProp = 2;

HaasenKinetics::HaasenKinetics(PolyCrystalMatT& poly) :
  SlipKinetics(poly)
{
  // fetch input file
  ifstreamT& in = poly.Input_x();

  // allocate space for material properties
  fMatProp.Dimension(kNumMatProp);

  // read material properties
  in >> fMatProp[0];     // "m" strain rate sensitivity exponent
  in >> fMatProp[1];     // tau_0
}

HaasenKinetics::~HaasenKinetics() { }

double HaasenKinetics::Phi(double tau, int is)
{
  // compute Phi
  double tausgn = fabs(tau)/tau;
  double qnt = ComputeInternalQnts(tau, is);
  return  qnt*tau*tausgn;
}

double HaasenKinetics::DPhiDTau(double tau, int is)
{
  // compute d(Phi)/d(Tau)
  double qnt = ComputeInternalQnts(tau, is);
  return  qnt/(fMatProp[0]*fMatProp[1]);
}

double HaasenKinetics::DPhiDIso(double tau, int is)
{
  // compute d(Phi)/d(Iso)
  double tausgn = fabs(tau)/tau;
  double qnt = ComputeInternalQnts(tau, is);
  return  -qnt/(fMatProp[0]*fMatProp[1])*tausgn;
}

double HaasenKinetics::DPhiDKin(double tau, int is)
{
  // compute d(Phi)/d(Kin)
  #pragma unused(tau, is)
  return 0.;
}

double HaasenKinetics::Psi(double gamdot, int is)
{
  // compute Psi  (not used)
  #pragma unused(gamdot, is)
  throwRunTimeError("HaasenKinetics::Psi called");
  return 0.;
}

double HaasenKinetics::DPsiDGamdot(double gamdot, int is)
{
  // compute d(Psi)/d(gamdot)
  #pragma unused(gamdot, is)
  throwRunTimeError("HaasenKinetics::DPsiDGamdot called");
  return 0.;
}

double HaasenKinetics::DPsiDIso(double gamdot, int is)
{
  // compute d(Psi)/d(Iso)
  #pragma unused(gamdot, is)
  throwRunTimeError("HaasenKinetics::DPsiDIso called");
  return 0.;
}

double HaasenKinetics::DPsiDKin(double gamdot, int is)
{
  // compute d(Psi)/d(Kin)
  #pragma unused(gamdot, is)
  throwRunTimeError("HaasenKinetics::DPsiDKin called");
  return 0.;
}

void HaasenKinetics::Print(ostream& out) const
{
  // print input values for kinetic equation of slip
  out << "       Rate sensitivity exponent (m) . . . . . . = " << fMatProp[0] << "\n";
  out << "       Tau_0 . . . . . . . . . . . . . . . . . . = " << fMatProp[1] << "\n";
}

void HaasenKinetics::PrintName(ostream& out) const
{
  // print model name
  out << "    Haasen's K.E. with nondirectional hardening\n";
}

double HaasenKinetics::ComputeInternalQnts(double& tau, const int is)
{
  double C = fHard.ComputeHardQnts(is);
  tau = (fabs(tau) - fHard.IsoHardeningStress(is))/fMatProp[1];
  if (tau <= 0.) tau = 0.;
  return  C * pow(tau, 1./fMatProp[0]-1.);
}
