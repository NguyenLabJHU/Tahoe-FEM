/*
  File: SlipKinetics.h

  General Form:  
               GammaDot = Phi(variables)

  General Inverse Form:
               Tau = Psi(variables)

  Currently:
      Phi = Gdot0*(|tau|/g)^(1/m)*sgn(tau)                Power Law I
      Phi = Gdot0*(|tau-alpha|/g)^(1/m)*sgn(tau-alpha)    Power Law II
      Phi = C*(|tau-alpha|-g)^(1/m)*sgn(tau-alpha)        Haasen Power Law
            with C = b*B(T)*MobileDislocDensity (alpha=0)
  where: g = iso, alpha = kin
*/

#ifndef _SLIP_KINETICS_H_
#define _SLIP_KINETICS_H_

#include "SlipHardening.h"
#include "dArrayT.h"


namespace Tahoe {

class PolyCrystalMatT;

class SlipKinetics
{
 public:
  // enum variable for crystal kinetic equation models
  enum SlipKinEqn { kPowLawI  = 1,     // dir hardening
                    kPowLawII = 2,     // dir/nondir hardening
                    kHaasen   = 3 };   // Haasen's model

  // constructor
  SlipKinetics(PolyCrystalMatT& poly);

  // destructor
  virtual ~SlipKinetics();

  // kinetic equation and its derivatives
  virtual double Phi(double tau, int is) = 0;
  virtual double DPhiDTau(double tau, int is) = 0;
  virtual double DPhiDIso(double tau, int is) = 0;
  virtual double DPhiDKin(double tau, int is) = 0;

  // inverse of kinetic equation and its derivatives
  virtual double Psi(double gamdot, int is) = 0;
  virtual double DPsiDGamdot(double gamdot, int is) = 0;
  virtual double DPsiDIso(double gamdot, int is) = 0;
  virtual double DPsiDKin(double gamdot, int is) = 0;

  // accesor to material constants of kinetic equation
  const dArrayT& MaterialProperties() const; 

  // print kinetic equation data and model name
  virtual void Print(ostream& out) const = 0;
  virtual void PrintName(ostream& out) const = 0;

  // continuation method using the rate sensitivity exponent
  virtual void SetUpRateSensitivity();
  virtual void ComputeRateSensitivity();
  virtual bool IsMaxRateSensitivity();
  virtual void RestoreRateSensitivity();

 protected:
  // compute max/min values of power law argument
  void MaxMinArgPowerLaw(const double &xm);

  // power
  double Power(const double &x, const double &y);

  // checking argument of power law
  void CheckArgumentRange(double &arg, const double &sign);

 protected:
  // reference to hardening object
  SlipHardening& fHard;

  // material properties array
  dArrayT fMatProp;

  // max & min value of power law argument
  double fArgMin;
  double fArgMax;
};

inline const dArrayT& SlipKinetics::MaterialProperties() const
{ return fMatProp; }

} // namespace Tahoe 
#endif  /* _SLIP_KINETICS_H_ */
