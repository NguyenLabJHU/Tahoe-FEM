/*
  File: HaasenHardening.h
*/

#ifndef _HAASEN_HARDENING_H_
#define _HAASEN_HARDENING_H_

#include <iostream.h>

#include "SlipHardening.h"
#include "dArrayT.h"
#include "dMatrixT.h"


namespace Tahoe {

class PolyCrystalMatT;
class HaasenKinetics;
class ifstreamT;

class HaasenHardening : public SlipHardening
{
 public:
  // constructor
  HaasenHardening(PolyCrystalMatT& poly);

  // destructor
  ~HaasenHardening();

  // initialize hardening variables at t_0
  virtual void InitializeHardVariables();

  // magnitude of hardening array
  virtual double Magnitude() const;

  // total number of hardening variables
  virtual const int NumberOfVariables() const;

  // update/reset history of hardening variables
  virtual void UpdateHistory();
  virtual void ResetHistory();

  // fetch hardening variables
  virtual void LoadHardData(int dim, int dex, dArrayT& d_array);

  // solution-related methods
  virtual void ExplicitUpdateHard();               // forward Euler estimate
  virtual void ImplicitUpdateHard();               // implicit estimate
  virtual void ImplicitSolveHard();                // backward Euler estimate
  virtual void FormRHS(const dArrayT& tauIso, dArrayT& rhs);     // residual
  virtual void FormLHS(const dArrayT& tauIso, dMatrixT& lhs);    // jacobian
  virtual bool Converged(double toler);            // check solution convergence

  // hardening modulus 
  virtual double HardeningModulus() const;

  // accesors to hardening quantities
  virtual const dArrayT& ComputeHardQnts();
  virtual const double ComputeHardQnts(int is);

  // print hardening data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 private:
  // nondirectional slip system hardness
  void NonDirectionalHardening(dArrayT& tauIso, const dArrayT& DDtot,
			       const dArrayT& Xaa);

  // temperature dependent coefficients
  void TempDependentFuncs(ifstreamT& in);

  // slip system interaction coefficients
  void SlipInteractionMatrix(ifstreamT& in);

  // latent hardening matrix
  void LatentHardMatrix();

  // compute some constant/variable hardening quantities
  void InternalHardQntsConst();
  void InternalHardQntsVaria(dArrayT& DDtot, dArrayT& Xaa);

  // hardening laws for DDtot, DDimm, and Xaa
  const double HardeningLawDDtot(double DDtot, int iss);
  const double HardeningLawDDimm(double DDimm, int iss);
  const double HardeningLawXaa(double Xaa, int iss);

 private:
  // temperature dependent constants
  double fBT;
  double fGT;
  double fWT;

  // modifiers of Bassani-Wu interaction coefficients
  double fC1xFab;
  double fC2xFab;

  // norms to check solution
  double fNormDDt0;
  double fNormDDi0;
  double fNormXaa0;
  double fNormDDt;
  double fNormDDi;
  double fNormXaa;

  // coefficients in evolution eqns for hardening variables
  dArrayT fAaa;        // Xaa's evolution eqn
  dArrayT fC0aa;       // Xaa's evolution eqn
  dArrayT fCaa;        // Xaa's evolution eqn
  dArrayT fDa;         // DDimm's evolution eqn

  // Haasen hardening variables at t and t_n
  dArrayT fDDtot;      // total dislocation density
  dArrayT fDDimm;      // immobile dislocation density
  dArrayT fXaa;        // rearrangement coefficients
  dArrayT fDDtot_n;
  dArrayT fDDimm_n;
  dArrayT fXaa_n;

  // total shear strain at t and t_n
  dArrayT fGamma;
  dArrayT fGamma_n;

  // effective stress
  dArrayT fTauEff;

  // slip system interaction matrix
  dMatrixT fFab;

  // latent hardening matrix
  dMatrixT fLatentMatx;

  // workspace
  dArrayT farray;
};

} // namespace Tahoe 
#endif /* _HAASEN_HARDENING_H_ */
