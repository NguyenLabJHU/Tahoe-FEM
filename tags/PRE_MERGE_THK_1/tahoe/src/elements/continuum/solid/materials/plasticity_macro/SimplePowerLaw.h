/*
  File: SimplePowerLaw.h

  Kinetic Equation Form: (s=eqp)

    eqpdot = f(sigma, eqp) = edot0*[(sigma/syield0)^(1/m) - 1];
    syield0 = g(eqp) = s0*(eqp/e0 + 1)^n

    material properties: (s0, e0, edot0, n, m)
*/

#ifndef _SIMPLE_POWER_LAW_H_
#define _SIMPLE_POWER_LAW_H_

#include "KineticEqnBase.h"
#include "ifstreamT.h"
#include "dArrayT.h"


namespace Tahoe {

class EVPFDBaseT;

class SimplePowerLaw : public KineticEqnBase
{
 public:

  // constructor
  SimplePowerLaw(EVPFDBaseT& model);

  // destructor
  ~SimplePowerLaw();

  // static yield stress
  virtual double g(double eqp);

  // kinetic equation functions
  virtual double f        (double sigma, double eqp);
  virtual double DfDsigma (double sigma, double eqp);
  virtual double DfDs     (double sigma, double eqp);

  // dynamic yield condition functions
  virtual double h         (double eqpdot, double kappa);
  virtual double DhDeqpdot (double eqpdot, double kappa);
  virtual double DhDs      (double eqpdot, double kappa);

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 private:

};

} // namespace Tahoe 
#endif  /* _SIMPLE_POWER_LAW */
