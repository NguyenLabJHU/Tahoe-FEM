/*
  File: BCJKineticEqn.h

  Kinetic Equation:

    eqpdot = f(s_eff, kappa) = f * sinh[(s_eff - kappa - Y)/V]

  Dynamic Yield Condition:

    s_eff = h(eqpdot, kappa) = kappa + Y + V * [sinh(eqpdot/f)]^(-1)

  Material parameters: f, Y, V

          V = C1*exp(-C2/Theta)
          Y = C3/(C21+exp(-C4/Theta))*0.5*(1+tanh(C19(C20-Theta)))
          f = C5*exp(-C6/Theta)

    with 
       C1,C2,C3,C4,C5,C6 : material constants
       C19,C20,C21       : material constants
       Theta : temperature
*/

#ifndef _BCJ_KINETIC_EQN_H_
#define _BCJ_KINETIC_EQN_H_

#include "KineticEqnBase.h"
#include "dArrayT.h"


namespace Tahoe {

class EVPFDBaseT;

class BCJKineticEqn : public KineticEqnBase
{
 public:
  // constructor  
  BCJKineticEqn(EVPFDBaseT& model);

  // destructor
  ~BCJKineticEqn();

  // static yield stress
  virtual double g(double eqp);

  // kinetic equation functions
  virtual double f        (double sigma, double kappa);
  virtual double DfDsigma (double sigma, double kappa);
  virtual double DfDs     (double sigma, double kappa);

  // dynamic yield condition functions
  virtual double h         (double eqpdot, double kappa);
  virtual double DhDeqpdot (double eqpdot, double kappa);
  virtual double DhDs      (double eqpdot, double kappa);

  // print data and model name
  virtual void Print     (ostream& out) const;
  virtual void PrintName (ostream& out) const;

 private:
  // material properties for Kinetic Eqn
  void ComputeMaterialProperties(const double theta);

 private:
  // temperature
  double fTheta;

  // material constants
  double fC1, fC2, fC3, fC4, fC5, fC6;
  double fC19, fC20, fC21;
};

} // namespace Tahoe 
#endif  /* _BCJ_KINETIC_EQN_H_ */
