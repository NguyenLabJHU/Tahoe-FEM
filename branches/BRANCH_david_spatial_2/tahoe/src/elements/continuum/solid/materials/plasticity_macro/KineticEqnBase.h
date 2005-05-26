/*
  File: KineticEqnBase.h
*/

#ifndef _KINETIC_EQN_BASE_H_
#define _KINETIC_EQN_BASE_H_

#include "dArrayT.h"


namespace Tahoe {

class KineticEqnBase
{
 public:

  // enum variable to select kinetic equation model
  enum KineticEqnModel { kSimplePowLaw = 1,
                         kBCJKinEqn    = 2 };  // BCJ Model only

  // constructor
  KineticEqnBase() {;}

  // destructor
  virtual ~KineticEqnBase() {;}

  // static yield stress
  virtual double g(double s) = 0;
  
  // kinetic equation functions
  virtual double f        (double sigma, double s) = 0;
  virtual double DfDsigma (double sigma, double s) = 0;
  virtual double DfDs     (double sigma, double s) = 0;

  // dynamic yield condition functions
  virtual double h         (double eqpdot, double s) = 0;
  virtual double DhDeqpdot (double eqpdot, double s) = 0;
  virtual double DhDs      (double eqpdot, double s) = 0;

 protected:
  // array of material properties
  dArrayT fMatProp;
};

//KineticEqnBase::KineticEqnBase() { }
//KineticEqnBase::~KineticEqnBase() { }

} // namespace Tahoe 
#endif  /* _KINETIC_EQN_BASE_H_ */
