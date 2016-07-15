/*
  File: SinhKinetics.cpp
*/

#include "SinhKinetics.h"
#include "PolyCrystalMatT.h"
#include "Utils.h"


using namespace Tahoe;

const int kNumMatProp = 2;

SinhKinetics::SinhKinetics(PolyCrystalMatT& poly) :
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

SinhKinetics::~SinhKinetics() { }

double SinhKinetics::Phi(double tau, int is)
{
  // compute Phi
  double tauIso = fHard.IsoHardeningStress(is);
  double tmp = tau/tauIso;

  double qnt = sinh(tmp);
  return  fMatProp[1]*qnt;
}

double SinhKinetics::DPhiDTau(double tau, int is)
{
  // compute d(Phi)/d(Tau)
  double tauIso = fHard.IsoHardeningStress(is);
  double tmp = tau/tauIso;

  double qnt = cosh(tmp);
  return  fMatProp[1]/tauIso*qnt;
}

double SinhKinetics::DPhiDIso(double tau, int is)
{
  // compute d(Phi)/d(Iso)
  double tauIso = fHard.IsoHardeningStress(is);
  double tmp = tau/tauIso;

  double qnt = cosh(tmp);
  return  -fMatProp[1]*tmp/(tauIso)*qnt;
}

double SinhKinetics::DPhiDKin(double tau, int is)
{
  // compute d(Phi)/d(Kin)
  #pragma unused(tau, is)
  return 0.;
}

double SinhKinetics::Psi(double gamdot, int is)
{
  // compute Psi
  double tauIso = fHard.IsoHardeningStress(is);
    double tmp = gamdot/fMatProp[1];
    double qnt = asinh(tmp);
  return  tauIso*qnt;
}

double SinhKinetics::DPsiDGamdot(double gamdot, int is)
{
  // compute d(Psi)/d(gamdot)
  double tauIso = fHard.IsoHardeningStress(is);
    double tmp = gamdot/fMatProp[1];
    double qnt = sqrt(1.0+tmp*tmp);
  return  (tauIso/fMatProp[1]/qnt);
}

double SinhKinetics::DPsiDIso(double gamdot, int is)
{
    double tauIso = fHard.IsoHardeningStress(is);
    double tmp = gamdot/fMatProp[1];
    double qnt = asinh(tmp);
    return  qnt;
}

double SinhKinetics::DPsiDKin(double gamdot, int is)
{
  // compute d(Psi)/d(Kin)
  #pragma unused(gamdot, is)
  return  0.;
}

