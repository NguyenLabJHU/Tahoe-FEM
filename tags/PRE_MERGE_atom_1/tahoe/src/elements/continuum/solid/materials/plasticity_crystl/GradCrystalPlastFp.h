/* $Id: GradCrystalPlastFp.h,v 1.4 2002-11-14 17:06:32 paklein Exp $ */
#ifndef _GRAD_CRYSTAL_PLAST_FP_H_
#define _GRAD_CRYSTAL_PLAST_FP_H_

#include "LocalCrystalPlastFp.h"
#include "CrystalElasticity.h"
#include "GradientTools.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

namespace Tahoe {

class GradCrystalPlastFp : public LocalCrystalPlastFp
{
 public:
  // constructor
  GradCrystalPlastFp(ifstreamT& in, const FDMatSupportT& support);

  // destructor
  ~GradCrystalPlastFp();

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

  // update/reset crystal state
  virtual void UpdateHistory();
  virtual void ResetHistory();

  // output related methods
  virtual int NumOutputVariables() const;
  virtual void OutputLabels(ArrayT<StringT>& labels) const;
  virtual void ComputeOutput(dArrayT& output);

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

 protected:
  // slip kinetics
  virtual void SetSlipKinetics();

  // slip hardening law
  virtual void SetSlipHardening();

  // recover crystal variables
  virtual void LoadCrystalData(ElementCardT& element, int intpt, int igrain);

  // solve for the state at each crystal
  virtual void SolveCrystalState();

  // solve for incremental shear strain
  virtual void SolveForPlasticDefGradient(int& ierr);

  // slip system resolve shear & back stresses
  virtual void ResolveShearStress();

 private:

  // number of crystal variables to be stored
  virtual int NumVariablesPerElement(void);

  // initial value of crystal variables
  virtual void InitializeCrystalVariables(ElementCardT& element);

  // fetch crystal curvature and associated stress
  void LoadCrystalCurvature(ElementCardT& element, int intpt, int igrain);

  // convergence of crystal state
  bool CheckConvergenceOfState() const;

  // dislocation tensor aBar in Bbar configuration
  void LatticeCurvature(ElementCardT& element, int igrn);

  // add gradient dependent term for Local Jacobian
  virtual void AddGradTermToLHS(dMatrixT& lhs, const dMatrixT& matx);

  // add gradient dependent term to consistent moduli
  virtual void AddGradTermToC_ijkl();

 protected:
  // number of element vertex nodes for spatial gradient evaluation
  // assumes: Quad in 2D (fNumNodes = 4); Hexa in 3D (fNumNodes = 8)
  int fNumNodes;

  // refs to nodal initial coords of element
  const LocalArrayT& fLocInitX;

  // nodal coords at current configuration
  LocalArrayT fLocCurrX;

  // pointer to supporting class for gradient evaluations
  GradientTools* fGradTool;

  // plastic deformation gradients at IPs and nodes
  ArrayT<dMatrixT> fFpIP;
  ArrayT<dMatrixT> fFpNodes;

  // spatial gradient of plastic deformation gradient
  ArrayT<dMatrixT> fGradFp;

  // the Curl of Fp^T
  dMatrixT fCurlFpT;

  // lattice curvatures and associated stress
  dMatrixT fKe_n;
  dMatrixT fKe;
  dMatrixT fXe;

  // workspaces for norms of Fp and hardness
  dArrayT fnormFp;
  dArrayT fnormFp0;
  dArrayT fnormHard;
  dArrayT fnormHard0;

  // worspace
  dMatrixT fMatx4;
};

} // namespace Tahoe 
#endif /* _GRAD_CRYSTAL_PLAST_FP_H_ */
