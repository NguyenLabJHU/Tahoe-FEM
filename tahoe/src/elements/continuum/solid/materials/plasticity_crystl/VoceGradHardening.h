/*
  File: VoceGradHardening.h
*/

#ifndef _VOCE_GRAD_HARDENING_H_
#define _VOCE_GRAD_HARDENING_H_

#include <iostream.h>

#include "SlipHardening.h"
#include "dArrayT.h"
#include "dMatrixT.h"


namespace Tahoe {

class PolyCrystalMatT;

class VoceGradHardening : public SlipHardening
{
 public:
  // constructor
  VoceGradHardening(PolyCrystalMatT& poly);

  // destructor
  ~VoceGradHardening();

  // initialize hardening variables at t_0
  virtual void InitializeHardVariables();

  // magnitude of hardness array
  virtual double Magnitude () const;

  // total number of hardening variables
  virtual const int NumberOfVariables() const;

  // update/reset history of hardening variables
  virtual void UpdateHistory();
  virtual void ResetHistory();

  // fetch hardening variables
  virtual void LoadHardData(int dim, int dex, dArrayT& d_array);

  // solution-related methods
  virtual void ExplicitUpdateHard();            // forward Euler estimate
  virtual void ImplicitUpdateHard();            // implicit estimate
  virtual void ImplicitSolveHard();             // backward Euler estimate
  virtual void FormRHS(const dArrayT& tauIso, dArrayT& rhs);  // residual
  virtual void FormLHS(const dArrayT& tauIso, dMatrixT& lhs); // Jacobian
  virtual bool Converged(double toler);         // check convergence

  // hardening modulus 
  virtual double HardeningModulus() const;

  // overides accesor to isotropic hardening variable in base class
  virtual const double IsoHardeningStress (int is) const;

 private:
  // compute some hardening quantities
  void InternalHardQnts();

  // hardening law
  const double HardeningLaw(double tauIso, int kcode);

 private:
  // saturation value of hardness
  double fTauIsoSat;

  // isotropic hardening stress at t_n (SS dislocations)
  dArrayT fTauIso_n;

  // kinematic hardening stress (GN dislocations)
  // dArrayT fTauKin;  // define in SlipHardening class

  // some internal variables (shear and work rate)
  dArrayT fInternal;
};

} // namespace Tahoe 
#endif /* _VOCE_GRAD_HARDENING_H_ */
