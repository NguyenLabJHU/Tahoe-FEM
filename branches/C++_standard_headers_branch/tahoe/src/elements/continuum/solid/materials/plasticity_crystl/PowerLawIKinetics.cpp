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

  // set Max/Min values of argument in power law
  MaxMinArgPowerLaw(fMatProp[0]);

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

//  double qnt = Power( fabs(tau)/tauIso, (1./fMatProp[0]-1.) );
//  return  fMatProp[1]*(tau/tauIso)*qnt;

  double tmp = tau/tauIso;
  double sign;
  if (fabs(tau) == 0.0) 
     sign = 1.0;
  else
     sign = fabs(tau)/tau;

  if ( !(fabs(tau) == 0.0) ) CheckArgumentRange(tmp, sign);

  double qnt = Power( fabs(tmp), (1./fMatProp[0]-1.) );
  return  fMatProp[1]*(tmp)*qnt;
}

double PowerLawIKinetics::DPhiDTau(double tau, int is)
{
  // compute d(Phi)/d(Tau)
  double tauIso = fHard.IsoHardeningStress(is);
//  double qnt = pow(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
//  double qnt = exp( (1./fMatProp[0]-1.)*log(fabs(tau)/tauIso) );

//  double qnt = Power( fabs(tau)/tauIso, (1./fMatProp[0]-1.) );
//  return  fMatProp[1]/(fMatProp[0]*tauIso)*qnt;

  double tmp = tau/tauIso;
  double sign;
  if (fabs(tau) == 0.0) 
     sign = 1.0;
  else
     sign = fabs(tau)/tau;

  if ( !(fabs(tau) == 0.0) ) CheckArgumentRange(tmp, sign);

  double qnt = Power( fabs(tmp), (1./fMatProp[0]-1.) );
  return  fMatProp[1]/(fMatProp[0]*tauIso)*qnt;
}

double PowerLawIKinetics::DPhiDIso(double tau, int is)
{
  // compute d(Phi)/d(Iso)
  double tauIso = fHard.IsoHardeningStress(is);
//  double qnt = pow(fabs(tau)/tauIso, 1./fMatProp[0]-1.);
//  double qnt = exp( (1./fMatProp[0]-1.)*log(fabs(tau)/tauIso) );

//  double qnt = Power( fabs(tau)/tauIso, (1./fMatProp[0]-1.) );
//  return  -fMatProp[1]/(fMatProp[0]*tauIso)*(tau/tauIso)*qnt;

  double tmp = tau/tauIso;
  double sign;
  if (fabs(tau) == 0.0) 
     sign = 1.0;
  else
     sign = fabs(tau)/tau;

  if ( !(fabs(tau) == 0.0) ) CheckArgumentRange(tmp, sign);

  double qnt = Power( fabs(tmp), (1./fMatProp[0]-1.) );
  return  -fMatProp[1]/(fMatProp[0]*tauIso)*(tmp)*qnt;
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

void PowerLawIKinetics::SetUpRateSensitivity()
{
  if (fkmax > 50.e0) 
     fk = 2.5e0;
  else
     fk = fkmax;
}

void PowerLawIKinetics::ComputeRateSensitivity()
{
  fk = min(2.e0*fk, fkmax);
  fMatProp[0] = 1.e0 / fk;
  if (fk == fkmax) fMatProp[0] = fxm;

  // compute new Max/Min values of argument in power law
  MaxMinArgPowerLaw(fMatProp[0]);
}

bool PowerLawIKinetics::IsMaxRateSensitivity()
{
  return (fk == fkmax);
}

void PowerLawIKinetics::RestoreRateSensitivity()
{
   fMatProp[0] = fxm;

  // restore Max/Min values of argument in power law
  MaxMinArgPowerLaw(fMatProp[0]);
}
