/*
  File: PowerLawIKinetics.cpp
*/

#include "PowerLawIKinetics.h"
#include "PolyCrystalMatT.h"
#include "Utils.h"


using namespace Tahoe;

const int kNumMatProp = 2;

PowerLawIKinetics::PowerLawIKinetics(PolyCrystalMatT& poly) :
  SlipKinetics(poly)
{
  // fetch input file
  ifstreamT& in = poly.Input_x();

  // allocate space for material properties
  fMatProp.Dimension(kNumMatProp);

  // read material properties
  in >> fMatProp[0];     // "m" strain rate sensitivity exponent
  in >> fMatProp[1];     // "Gdot_0"

  // set up parameters for continuation method using "m"
  fxm   = fMatProp[0];
  fkmax = 1.e0 / fxm;
}

PowerLawIKinetics::~PowerLawIKinetics() { }

double PowerLawIKinetics::Phi(double tau, int is)
{
  // compute Phi
  double tauIso = fHard.IsoHardeningStress(is);
//  double qnt = pow(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
//  double qnt = exp( (1./fMatProp[0]-1.)*log(fabs(tau)/tauIso) );
  double qnt = Power( fabs(tau)/tauIso, (1./fMatProp[0]-1.) );
  return  fMatProp[1]*(tau/tauIso)*qnt;
}

double PowerLawIKinetics::DPhiDTau(double tau, int is)
{
  // compute d(Phi)/d(Tau)
  double tauIso = fHard.IsoHardeningStress(is);
//  double qnt = pow(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
//  double qnt = exp( (1./fMatProp[0]-1.)*log(fabs(tau)/tauIso) );
  double qnt = Power( fabs(tau)/tauIso, (1./fMatProp[0]-1.) );
  return  fMatProp[1]/(fMatProp[0]*tauIso)*qnt;
}

double PowerLawIKinetics::DPhiDIso(double tau, int is)
{
  // compute d(Phi)/d(Iso)
  double tauIso = fHard.IsoHardeningStress(is);
//  double qnt = pow(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
//  double qnt = exp( (1./fMatProp[0]-1.)*log(fabs(tau)/tauIso) );
  double qnt = Power( fabs(tau)/tauIso, (1./fMatProp[0]-1.) );
  return  -fMatProp[1]/(fMatProp[0]*tauIso)*(tau/tauIso)*qnt;
}

double PowerLawIKinetics::DPhiDKin(double tau, int is)
{
  // compute d(Phi)/d(Kin)
  #pragma unused(tau, is)
  return 0.;
}

double PowerLawIKinetics::Psi(double gamdot, int is)
{
  // compute Psi
  double tauIso = fHard.IsoHardeningStress(is);
  //double qnt = pow(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  double qnt = Power( fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  return  tauIso*(gamdot/fMatProp[1])*qnt;
}

double PowerLawIKinetics::DPsiDGamdot(double gamdot, int is)
{
  // compute d(Psi)/d(gamdot)
  double tauIso = fHard.IsoHardeningStress(is);
  //double qnt = pow(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  double qnt = Power( fabs(gamdot)/fMatProp[1], fMatProp[0]-1. );
  return  (fMatProp[0]*tauIso)/fMatProp[1]*qnt;
}

double PowerLawIKinetics::DPsiDIso(double gamdot, int is)
{
  // compute d(Psi)/d(Iso)
  #pragma unused(is)
  //double qnt = pow(fabs(gamdot)/fMatProp[1], fMatProp[0]-1.);
  double qnt = Power( fabs(gamdot)/fMatProp[1], fMatProp[0]-1. );
  return  gamdot/fMatProp[1]*qnt;
}

double PowerLawIKinetics::DPsiDKin(double gamdot, int is)
{
  // compute d(Psi)/d(Kin)
  #pragma unused(gamdot, is)
  return  0.;
}

void PowerLawIKinetics::Print(ostream& out) const
{
  // print input values for kinetic equation of slip
  out << "       Rate sensitivity exponent (m) . . . . . . = " << fMatProp[0] << "\n";
  out << "       Gdot_0. . . . . . . . . . . . . . . . . . = " << fMatProp[1] << "\n";
}

void PowerLawIKinetics::PrintName(ostream& out) const
{
  // print model name
  out << "    Power law K.E. with nondirectional hardening\n";
}

void PowerLawIKinetics::SetUpRateSensitivity()
{
  if (fkmax > 30.e0) 
     fk = 2.5e0;
  else
     fk = fkmax;
}

void PowerLawIKinetics::ComputeRateSensitivity()
{
  fk = min(2.e0*fk, fkmax);
  fMatProp[0] = 1.e0 / fk;
  if (fk == fkmax) fMatProp[0] = fxm;
}

bool PowerLawIKinetics::IsMaxRateSensitivity()
{
   return (fk == fkmax);
}

void PowerLawIKinetics::RestoreRateSensitivity()
{
   fMatProp[0] = fxm;
}
