/*
  File: BCJHypo3D.h
*/

#ifndef _BCJ_HYPO_3D_H_
#define _BCJ_HYPO_3D_H_

#include "EVPFDBaseT.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "SpectralDecompT.h"

class ifstreamT;
class ElasticT;
class ElementCardT;
class StringT;

class BCJHypo3D : public EVPFDBaseT
{
 public:
  // constructor
  BCJHypo3D(ifstreamT& in, const FiniteStrainT& element);

  // destructor
  ~BCJHypo3D();

  // number of variables to be stored
  virtual int NumVariablesPerElement();

  // Cauchy stress
  virtual const dSymMatrixT& s_ij();   

  // tangent modulus
  virtual const dMatrixT& c_ijkl();

  // form residual
  virtual void FormRHS(const dArrayT& variab, dArrayT& rhs);

  // form Jacobian
  virtual void FormLHS(const dArrayT& variab, dMatrixT& lhs);

  // update/reset state
  virtual void UpdateHistory();
  virtual void ResetHistory();

  // output related methods
  virtual int NumOutputVariables() const;
  virtual void OutputLabels(ArrayT<StringT>& labels) const;
  virtual void ComputeOutput(dArrayT& output);

  // print data and model name
  virtual void Print(ostream& out) const;
  virtual void PrintName(ostream& out) const;

  /* form of tangent matrix */
  virtual GlobalT::SystemTypeT TangentType(void) const;

 protected:
  // indexes to access internal variable (scalars) array
  enum InternalVariables { kDEQP = 0,
			   kALPH = 1,
			   kKAPP = 2 };

  enum EQValues { kEQP_n   = 0,
		  kEQP     = 1,
		  kEQSig_n = 2,
		  kEQSig   = 3,
		  kPresTr  = 4 };

  // material properties for isotropic/kinematic hardening rules
  void ComputeMaterialProperties(double theta);

  // size of system of nonlinear equations
  virtual int GetNumberOfEqns();

  // kinetic equation
  virtual void SetKineticEquation();

  // initial value of variables
  virtual void InitializeVariables();

  // recover variables
  virtual void LoadElementData(ElementCardT& element, int intpt);

  // solve for the state at each IP
  virtual void SolveState();

  // backward integration of constitutive equations
  virtual void IntegrateConstitutiveEqns(bool& converged, int subIncr);

  // polar decomposition of relative deformation gradient
  void PolarDecomposition();

  // incremental total strain (logaritmic strain)
  virtual void IncrementalStrain();

  // rotate tensorial state variables at t_n to current configuration
  virtual void RotateTensorialVariables();

  // trial stresses
  virtual void ElasticTrialStress();

  // solve for state variables
  virtual void Solve();

  // Cauchy stress and backstress
  virtual void UpdateStresses();

  // tangent moduli
  virtual void TangentModuli();

  // elastic moduli
  virtual void ElasticModuli();

 private:
  // forward gradient estimate for DEQP, DALP, DKAPP
  void ForwardGradientEstimate();

  // internal quantities
  void ComputeInternalQntsRHS(const dArrayT& array);
  void ComputeInternalQntsLHS(const dArrayT& array);

 protected:
  // material constants
  double fC7, fC8, fC9, fC10, fC11, fC12;
  double fC13, fC14, fC15, fC16, fC17, fC18;

  // scalar stress variables
  double fEQSigTr;
  double fradial;     // fradial=fEQValues[EQSig]/fEQSigTr
  double ffactor;

  // other scalar variables
  double fEta;
  double fDEtaDDEQP;
  double fDEtaDALPH;
  double fXMag;

  // material properties
  dArrayT fMatProp;

  // deformation gradients
  dMatrixT fFtot_n;   // at t_n
  dMatrixT fF;        // relative def gradient = Ftot*Ftot_n^(-1)
  dMatrixT fFr;       // relative def gradient in continuation method

  // backstress
  dSymMatrixT falph_ij;

  // Cauchy stress and backstress at t_n
  dSymMatrixT fs_ij_n;
  dSymMatrixT falph_ij_n;

  // Cauchy stress and backstress at t_n rotated to current configuration
  dSymMatrixT fsigma_n;
  dSymMatrixT falpha_n;

  // trial stresses
  dSymMatrixT fSigTr;
  dSymMatrixT fSigTrDev;
  dSymMatrixT fXiTr;

  // incremental strains
  dSymMatrixT fDEBar;
  dSymMatrixT fDE;

  // helpful tensors
  dSymMatrixT fA;
  dSymMatrixT fX;
  dSymMatrixT fUnitNorm;
  dSymMatrixT fUnitM;

  // second order identity tensors
  dSymMatrixT fISym;
  dMatrixT fI;

  // spectral/polar decomposition of fF
  SpectralDecompT fSpecD;
  dArrayT fEigs;
  dSymMatrixT fC;
  dSymMatrixT fU;
  dMatrixT fR;

  // array for scalar internal variables
  dArrayT fInternal_n;    // DEQP_n, ALPH_n, KAPP_n
  dArrayT fInternal;      // DEQP, ALPH, KAPP
  dArrayT fInt_save;
  dArrayT fEQValues;      // EQP_n, EQP, EQSig_n, EQSig, PresTr

  // arrays used in forward gradient estimate
  dArrayT fRHS;
  dMatrixT fLHS;

  // general workspaces
  dArrayT farray;
  dMatrixT fmatx1;
  dMatrixT fmatx2;
  dMatrixT fmatx3;
  dMatrixT fmatx4;
  dMatrixT fRank4;
  dMatrixT fIdentity4;
  dSymMatrixT fsymmatx1;
};

#endif /* _BCJ_HYPO_3D_ */
