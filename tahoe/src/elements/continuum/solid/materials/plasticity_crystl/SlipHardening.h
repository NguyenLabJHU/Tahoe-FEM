/*
  File: SlipHardening.h
*/

#ifndef _SLIP_HARDENING_H_
#define _SLIP_HARDENING_H_

#include <iostream.h>

#include "NLCSolverWrapperPtr.h"
#include "dArrayT.h"


namespace Tahoe {

class PolyCrystalMatT;
class GradCrystalPlast;
class GradCrystalPlastFp;
class PowerLawIKinetics;
class PowerLawIIKinetics;
class HaasenKinetics;
class NLCSolver;
class ifstreamT;
class dMatrixT;

class SlipHardening
{
  friend class GradCrystalPlast;
  friend class GradCrystalPlastFp;

 public:
  // enumeration for hardening models for local crystal plast
  enum LocHardModel { kHard_L1 = 1,     // iso hardening type (Voce's model)  
		      kHard_L2 = 2,     // iso/kin hardening type
                      kHard_L3 = 3 };   // iso hardening type (Haasen's model)

  // enumeration for hardening models for gradient crystal plast
  enum GradHardModel { kHard_G1 = 11 }; // ss-gn dislocation density type 

  // constructor
  SlipHardening(PolyCrystalMatT& poly);

  // destructor
  virtual ~SlipHardening();

  // nonlinear solver for hardening model
  void SetHardeningSolver(ifstreamT& in, int numvar);

  // set hardening variables at t_0
  virtual void InitializeHardVariables() = 0;

  // magnitude of hardening array
  virtual double Magnitude() const = 0;

  // accesors
  virtual const int NumberOfVariables() const = 0;
  const dArrayT& MaterialProperties() const;

  // update/reset history of hardening variables
  virtual void UpdateHistory() = 0;
  virtual void ResetHistory() = 0;

  // fetch hardening variables
  virtual void LoadHardData(int dim, int dex, dArrayT& d_array) = 0;

  // solution-related methods
  virtual void ExplicitUpdateHard() = 0;               // forward Euler estimate
  virtual void ImplicitUpdateHard() = 0;               // implicit estimate
  virtual void ImplicitSolveHard() = 0;                // backward Euler estimate
  virtual void FormRHS(const dArrayT& variab, dArrayT& rhs) = 0;  // residual
  virtual void FormLHS(const dArrayT& variab, dMatrixT& lhs) = 0; // jacobian
  virtual bool Converged(double toler) = 0;           // check convergence

  virtual void SaveCurrentSolution() { };
  virtual void RestoreSavedSolution() { };

  // hardening modulus 
  virtual double HardeningModulus() const = 0;

  // compute needed hardening quantities
  virtual const dArrayT& ComputeHardQnts();
  virtual const double ComputeHardQnts(int is);

  // accessors to isotropic and kinematic hardening stresses
  const dArrayT& IsoHardeningStress () const;
  const dArrayT& KinHardeningStress () const;
  virtual const double IsoHardeningStress (int is) const;
  virtual const double KinHardeningStress (int is) const;

 protected:
  // reference to PolyCrystalMat object
  PolyCrystalMatT& fPolyXtal;

  // time step
  const double& fdt;

  // number of hardening variables
  int fNumHardVar;

  // resolved shear stress and incremental shear strain
  const dArrayT& fTau;
  const dArrayT& fDGamma;

  // isotropic(nondirectional)/kinematic(directional) hardening stresses at t
  dArrayT fTauIso;
  dArrayT fTauKin;

  // incompatibility stress-like quantity (nonlocal model)
  dArrayT fTauInc;

  // slip hardening material properties
  dArrayT fMatProp;

  // initial values of hardening variables
  dArrayT fInitHardValues;

  // pointer and handler to nonlinear solver
  NLCSolver* fSolver;
  NLCSolverWrapperPtr fSolverPtr;
};

inline const dArrayT& SlipHardening::MaterialProperties() const { return fMatProp; }
inline const dArrayT& SlipHardening::IsoHardeningStress() const { return fTauIso; }
inline const dArrayT& SlipHardening::KinHardeningStress() const { return fTauKin; }

} // namespace Tahoe 
#endif /* _SLIP_HARDENING_H_ */

