/*
  File: PowerLawIIKinetics.h
*/

#ifndef _POWER_LAW_II_KINETICS_H_
#define _POWER_LAW_II_KINETICS_H_

#include "SlipKinetics.h"


namespace Tahoe {

class PolyCrystalMatT;

class PowerLawIIKinetics: public SlipKinetics
{
 public:
  // constructor
  PowerLawIIKinetics(PolyCrystalMatT& poly);

  // destructor
  ~PowerLawIIKinetics();

  // power law equation and its derivatives
  virtual double Phi(double tau, int is);
  virtual double DPhiDTau(double tau, int is);
  virtual double DPhiDIso(double tau, int is);
  virtual double DPhiDKin(double tau, int is);

  // inverse of power law equation and its derivatives
  virtual double Psi(double gamdot, int is);
  virtual double DPsiDGamdot(double gamdot, int is);
  virtual double DPsiDIso(double gamdot, int is);
  virtual double DPsiDKin(double gamdot, int is);

  // print kinetic equation data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 private:
  // compute some preliminary quantities
  double ComputeInternalQnts(double& tau, double tauIso, int is);
};

} // namespace Tahoe 
#endif  /* _POWER_LAW_II_KINETICS_H_ */
