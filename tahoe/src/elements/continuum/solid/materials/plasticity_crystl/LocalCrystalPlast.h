/* $Id: LocalCrystalPlast.h,v 1.2 2001-07-03 01:35:35 paklein Exp $ */
/*
  File: LocalCrystalPlast.h
*/

#ifndef _LOCAL_CRYSTAL_PLAST_H_
#define _LOCAL_CRYSTAL_PLAST_H_

#include "PolyCrystalMatT.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "LAdMatrixT.h"

class ifstreamT;
class ElasticT;
class ElementCardT;
class StringT;

class LocalCrystalPlast : public PolyCrystalMatT
{
 public:
  // constructor
  LocalCrystalPlast(ifstreamT& in, const FiniteStrainT& element);

  // destructor
  ~LocalCrystalPlast();

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

  // form of tangent matrix
  virtual GlobalT::SystemTypeT TangentType(void) const;

  // some needed accesors (by slip hardening classes mainly)
  virtual const dArrayT& GetResolvedShearStress() const;
  virtual const dArrayT& GetIncrSlipShearStrain() const;

 protected:
  // slip kinetics
  virtual void SetSlipKinetics();

  // slip hardening law
  virtual void SetSlipHardening();

  // initial value of crystal variables
  virtual void InitializeCrystalVariables();

  // recover crystal variables
  virtual void LoadCrystalData(ElementCardT& element, int intpt, int igrain);

  // recover averaged stress and moduli at IP
  virtual void LoadAggregateData(ElementCardT& element, int intpt);

  // solve for the state at each crystal
  virtual void SolveCrystalState();
  virtual void IterateOnCrystalState(bool& stateConverged, int subIncr);

  // trial deformation quantities
  void TrialDeformation();

  // forward gradient approximation for dgamma
  void ForwardGradientEstimate();

  // check convergence on DGamma
  bool Converged(double toler);

  // crystal Cauchy stress
  virtual void CrystalS_ij();

  // elastic contribution to crystal moduli
  virtual void CrystalC_ijkl_Elastic();

  //plastic contribution to crystal moduli
  virtual void CrystalC_ijkl_Plastic();

  // solve for incremental shear strain
  virtual void SolveForDGamma() { };
  virtual void SolveForDGamma(int& ierr);

  // slip system resolve shear stress
  void ResolveShearStress(const dSymMatrixT& Ce);

  // inverse of incremental plastic deformation gradient
  virtual void DeltaFPInverse(const dArrayT& dgamma);

  // Second Piola Kirchhoff stress
  void CrystalPKIIStress(const dSymMatrixT& Ce);

  // fDFpi/dDGamma term in local Jacobian
  void dDFpidDGamma(const dArrayT& dgamma, ArrayT<dMatrixT>& array);

  // frequent term used to compute the crystal moduli
  void dTaudCe(const dMatrixT& Z, const dSymMatrixT& P, dSymMatrixT& symmatx);

  // update plastic deformation gradient
  void FPInverse();

  // polar decomposition of deformation gradient
  void PolarDecomp();

	// function to compute 3D deformations regardless of dimensionality of the
	// problem. For 2D, the out-of-plane direction is x3 and the deformation
	// is assumed to be plane strain
	void Compute_Ftot_3D(dMatrixT& F_3D) const;	
	void Compute_Ftot_3D(dMatrixT& F_3D, int ip) const;	
	void Compute_Ftot_last_3D(dMatrixT& F_3D) const;	

  // 4th order tensor: c_ijkl=0.5*(b_ik b_jl + b_il b_jk)
  void Set_I_b_Tensor(const dSymMatrixT& b, dMatrixT& c);

  // 4th order tensor transformation: Co_ijkl = F_iI F_jJ F_kK f_lL Ci_IJKL
  void FFFFC_3D(dMatrixT& Co, dMatrixT& Ci, const dMatrixT& F);

 protected:
  // number of hardening variables
  int fNumHard;

  // norms to check convergence on DGamma
  double fMagDGam0;
  double fMagDGam;

  // deformation gradient (continuation method)
  dMatrixT fFt;   

  // elastic deformation gradients
  dMatrixT fFeTr;
  dMatrixT fFe;

  // plastic deformation gradients
  dMatrixT fFpi_n;
  dMatrixT fFpi;
  dMatrixT fDFp;
  dMatrixT fDFpi;

  // tensors in intermediate configuration
  dSymMatrixT fCeBarTr;
  dSymMatrixT fCeBar;
  dSymMatrixT fEeBar;
  dSymMatrixT fSBar;

  // rotation matrix from Euler angles
  dMatrixT fRotMat;

  // crystal Cauchy stress
  dSymMatrixT fs_ij;
  
  // crystal (spatial) consistent tangent operator
  dMatrixT fc_ijkl;

  // anisotropic contribution to crystal elasticity matrix
  dMatrixT fCanisoLat;    // lattice frame
  dMatrixT fCanisoBar;    // intermediate configuration

  // Schmidt tensors in sample coords
  ArrayT<dMatrixT> fZ;
  ArrayT<dSymMatrixT> fP;

  // incremental slip system shearing rate
  dArrayT fDGamma_n;
  dArrayT fDGamma;
  dArrayT fdgam_save;

  // resolve shear stress on slip systems
  dArrayT fTau;

  // second order identity tensor
  dSymMatrixT fISym;

  // work spaces used in polar decomposition of fFe
  dArrayT fEigs;
  dMatrixT fRe;
  dSymMatrixT fUe;

  // general workspaces
  double fdeltaI;
  dMatrixT fmatx1;
  dMatrixT fmatx2;
  dMatrixT fmatx3;
  dMatrixT fmatx4;
  dMatrixT fRank4;
  dSymMatrixT fsymmatx1;
  dSymMatrixT fsymmatx2;

  // workspaces for computing moduli 
  LAdMatrixT fLHS;
  dArrayT fRHS;
  ArrayT<dSymMatrixT> fA;
  ArrayT<dSymMatrixT> fB;
  ArrayT<dMatrixT> farray;

  // mid-point rule parameters
  double fTheta;

  // equivalent stress
  dSymMatrixT fAvgStress;
};

#endif /* _LOCAL_CRYSTAL_PLAST_H_ */
