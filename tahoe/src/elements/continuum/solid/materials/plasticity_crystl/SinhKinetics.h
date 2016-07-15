/*
  File: SinhKinetics.h
*/

#ifndef _SINH_KINETICS_H_
#define _SINH_KINETICS_H_

#include "SlipKinetics.h"


namespace Tahoe {

class PolyCrystalMatT;

class SinhKinetics: public SlipKinetics
{
 public:
  // constructor
  SinhKinetics(PolyCrystalMatT& poly);

  // destructor
  ~SinhKinetics();

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


 private:

  // variables used in setting up continuation method
  double fxm;
  double fk;
  double fkmax;
};

} // namespace Tahoe 
#endif  /* _SINH_KINETICS_H_ */
