/*
  File: GradCrystalPlast.h
*/

#ifndef _GRAD_CRYSTAL_PLAST_H_
#define _GRAD_CRYSTAL_PLAST_H_

#include "LocalCrystalPlast.h"

#include "ArrayT.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"

class GradCrystalPlast : public LocalCrystalPlast
{
 public:
  // constructor
  GradCrystalPlast(ifstreamT& in, const ElasticT& element);

  // destructor
  ~GradCrystalPlast();

  // initialize arrays
  void Initialize();

  // number of crystal variables to be stored
  virtual int NumVariablesPerElement();

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

  // form residual
  virtual void FormRHS(const dArrayT& dgamma, dArrayT& rhs);

  // form Jacobian
  virtual void FormLHS(const dArrayT& dgamma, dMatrixT& lhs);

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

  // initial value of crystal variables
  virtual void InitializeCrystalVariables();

  // recover crystal variables
  virtual void LoadCrystalData(ElementCardT& element, int intpt, int igrain);

  // solve for the state at each crystal
  virtual void SolveCrystalState();

  // solve for incremental shear strain
  virtual void SolveForDGamma();

  // inverse of incremental plastic deformation gradient
  virtual void DeltaFPInverse(const dArrayT& dgamma);

//DEV - deprecated
  // deformation gradient
  //virtual const dMatrixT& DeformationGradient(const LocalArrayT& disp);

 private:
  // fetch crystal curvature and associated stress
  void LoadCrystalCurvature(ElementCardT& element, int intpt, int igrain);

  // derivatives of shape functions (to compute GradFe)
  void ShapeFunctionDeriv();

  // convergence of crystal state
  bool CheckConvergenceOfState() const;

  // dislocation tensor A2e in Bbar configuration
  void LatticeCurvature(ElementCardT& element, int igrn);

  // dFe/dx at integration points
  void ComputeGradFe(const dArray2DT& GDNa, const ArrayT<dMatrixT>& nodal, ArrayT<dMatrixT>& grad);

  void dKedDGamma(ElementCardT& element);
  virtual void CrystalC_ijkl_Plastic();
  void AddToLHS(dMatrixT& lhs);
  double HardFuncDerivative(double& dgamma, double& taus, double& taux, int kcode);

 protected:
  // counter for integration points (needed when computing F)
  int fIP;

  // number of element vertex nodes for spatial gradient evaluation
  // assumes: Quad in 2D (fNumNodes = 4); Hexa in 3D (fNumNodes = 8)
  int fNumNodes;

  // refs to nodal initial coords of element
  const LocalArrayT& fLocInitX;

  // nodal coords at current configuration
  LocalArrayT fLocCurrX;

  // elastic deformation gradients at IPs and nodes
  ArrayT<dMatrixT> fFeTrIP;
  ArrayT<dMatrixT> fFeIP;
  ArrayT<dMatrixT> fFeTrNodes;
  ArrayT<dMatrixT> fFeNodes;

  // arrays of shape function and their derivatives
  dArray2DT fLNa;
  ArrayT<dArray2DT> fLDNa;    // local derivatives (parental domain)
  ArrayT<dArray2DT> fGDNa;    // global derivatives (physical domain)

  // matrix for nodal extrapolation (least square smoothing)
  dArray2DT fNodalExtrap;

  // spatial gradients of elastic deformation gradient
  ArrayT<dMatrixT> fGradFeTr;
  ArrayT<dMatrixT> fGradFe;

  // displacement gradient
  dMatrixT fGradU;

  // lattice curvatures and associated stress
  dMatrixT fKe_n;
  dMatrixT fKe;
  dMatrixT fXe;

  // array to hold d(Ke)/d(DGamma)
  ArrayT<ArrayT<dMatrixT> > fdKe;

  // workspaces for norms of dgamma and hardness
  dArrayT fnormDGam0;
  dArrayT fnormHard0;
  dArrayT fnormDGam;
  dArrayT fnormHard;
};

#endif /* _GRAD_CRYSTAL_PLAST_H_ */
