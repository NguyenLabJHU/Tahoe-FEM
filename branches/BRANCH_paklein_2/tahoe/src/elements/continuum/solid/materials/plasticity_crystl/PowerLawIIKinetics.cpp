/*
  File: PowerLawIIKinetics.cpp
*/

#include "PowerLawIIKinetics.h"
#include "PolyCrystalMatT.h"


using namespace Tahoe;

const int kNumMatProp = 2;

PowerLawIIKinetics::PowerLawIIKinetics(PolyCrystalMatT& poly) :
  SlipKinetics(poly)
{
  // fetch input file
  ifstreamT& in = poly.Input_x();

  // allocate space for material properties
  fMatProp.Dimension(kNumMatProp);

  // read material properties
  in >> fMatProp[0];     // "m" strain rate sensitivity exponent
  in >> fMatProp[1];     // "Gdot_0"
}

PowerLawIIKinetics::~PowerLawIIKinetics() { }

double PowerLawIIKinetics::Phi(double tau, int is)
{
  // compute Phi
  double tauIso = fHard.IsoHardeningStress(is);
  double qnt = ComputeInternalQnts(tau, tauIso, is);
  return  fMatProp[1]*(tau/tauIso)*qnt;
}

double PowerLawIIKinetics::DPhiDTau(double tau, int is)
{
  // compute d(Phi)/d(Tau)
  double tauIso = fHard.IsoHardeningStress(is);
  double qnt = ComputeInternalQnts(tau, tauIso, is);
  return  fMatProp[1]/(fMatProp[0]*tauIso)*qnt;
}

double PowerLawIIKinetics::DPhiDIso(double tau, int is)
{
  // compute d(Phi)/d(Iso)
  double tauIso = fHard.IsoHardeningStress(is);
  double qnt = ComputeInternalQnts(tau, tauIso, is);
  return -fMatProp[1]/(fMatProp[0]*tauIso)*(tau/tauIso)*qnt;
}

double PowerLawIIKinetics::DPhiDKin(double tau, int is)
{
  // compute d(Phi)/d(Kin)
  double tauIso = fHard.IsoHardeningStress(is);
  double qnt = ComputeInternalQnts(tau, tauIso, is);
  return  -fMatProp[1]/(fMatProp[0]*tauIso)*qnt;
}

double PowerLawIIKinetics::Psi(double gamdot, int is)
{
  // compute Psi
  double tauIso = fHard.IsoHardeningStress(is);
  double tauKin = fHard.KinHardeningStress(is);
//  double qnt = pow(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  double qnt = Power(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  return  tauKin + tauIso*(gamdot/fMatProp[1])*qnt;
}

double PowerLawIIKinetics::DPsiDGamdot(double gamdot, int is)
{
  // compute d(Psi)/d(gamdot)
  double tauIso = fHard.IsoHardeningStress(is);
//  double qnt = pow(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  double qnt = Power(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  return  (fMatProp[0]*tauIso)/fMatProp[1] * qnt;
}

double PowerLawIIKinetics::DPsiDIso(double gamdot, int is)
{
  // compute d(Psi)/d(Iso)
  #pragma unused(is)
//  double qnt = pow(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  double qnt = Power(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  return  gamdot/fMatProp[1]*qnt;
}

double PowerLawIIKinetics::DPsiDKin(double gamdot, int is)
{
  // compute d(Psi)/d(Kin)
  #pragma unused(gamdot, is)
  return  1.;
}

void PowerLawIIKinetics::Print(ostream& out) const
{
  // print input values for kinetic equation of slip
  out << "       Rate sensitivity exponent (m) . . . . . = " << fMatProp[0] << "\n";
  out << "       Gdot_0. . . . . . . . . . . . . . . . . = " << fMatProp[1] << "\n";
}

void PowerLawIIKinetics::PrintName(ostream& out) const
{
  // print model name
  out << "    Power law K.E. with directional/nodirectional hardening\n";
}

double PowerLawIIKinetics::ComputeInternalQnts(double& tau, double tauIso, int is)
{
  tau = tau - fHard.KinHardeningStress(is);
//  return  pow(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
  return  Power(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
}
