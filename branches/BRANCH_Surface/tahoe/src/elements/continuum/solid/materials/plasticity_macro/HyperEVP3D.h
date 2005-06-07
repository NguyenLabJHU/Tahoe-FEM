/* $Id: HyperEVP3D.h,v 1.7 2004-07-15 08:29:14 paklein Exp $ */
#ifndef _HYPER_EVP_3D_H_
#define _HYPER_EVP_3D_H_

#include "EVPFDBaseT.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"
#include "SpectralDecompT.h"


namespace Tahoe {

class ifstreamT;
class SolidElementT;
class ElementCardT;
class StringT;

class HyperEVP3D : public EVPFDBaseT
{
 public:
  // constructor
  HyperEVP3D(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~HyperEVP3D();

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

  /* form of tangent matrix */
  virtual GlobalT::SystemTypeT TangentType(void) const;

 protected:
  // indexes to access internal variable array
  enum InternalVariables { kEQP_n   = 0,
			   kEQP     = 1,
			   kEQSig_n = 2,
			   kEQSig   = 3,
			   kpressTr = 4 };

  // size of system of nonlinear equations
  virtual int GetNumberOfEqns();

  // kinetic equation
  virtual void SetKineticEquation();

  // recover variables
  virtual void LoadElementData(ElementCardT& element, int intpt);

  // solve for the state at each integration point
  virtual void IntegrateConstitutiveEqns();

  // polar decomposition of deformation gradient
  void PolarDecomposition();

  // trial strain (logarithmic strain)
  virtual void LogarithmicStrain();

  // trial stresses
  virtual void ElasticTrialStress();

  // solve for state variables
  virtual void Solve();

  // Cauchy stress
  virtual void CauchyStress();

  // update plastic deformation gradient
  virtual void FPInverse();

  // tangent moduli
  virtual void TangentModuli();

  // elastic moduli
  virtual void ElasticModuli();

 private:

  // number of variables to be stored
  virtual int NumVariablesPerElement();

  // initial value of variables
  virtual void InitializeVariables(ElementCardT& element);

  // forward gradient estimate for DEQP
  void ForwardGradientEstimate();

 protected:
  // radial return factor
  double fEta;

  // trial equivalent stress
  double fEQSigTr;
	
  // elastic deformation gradients
  dMatrixT fFeTr;
  dMatrixT fFe;

  // plastic deformation gradients
  dMatrixT fFpi_n;
  dMatrixT fFpi;
  dMatrixT fDFpi;

  // tensors in intermediate configuration
  dSymMatrixT fCeBarTr;
  dSymMatrixT fEeBarTr;
  dSymMatrixT fSigBarTr;
  dSymMatrixT fSigBarTrDev;
  dSymMatrixT fSigBar;

  // second order identity tensor
  dSymMatrixT fISym;

  // spectral/polar decomposition
  SpectralDecompT fSpecD;
  dArrayT fEigs;
  dMatrixT fReTr;
  dSymMatrixT fUeTr;

  // array for scalar internal variables
  // (EQP_n, EQP, EQSig_n, EQSig, pressTr)
  dArrayT fInternal;

  // incremental equiv plastic strain
  dArrayT fDEQP;

  // general workspaces
  dArrayT farray;
  dMatrixT fmatx1;
  dMatrixT fRank4;
  dMatrixT fIdentity4;
  dSymMatrixT fsymmatx1;
  dSymMatrixT fsymmatx2;
};

} // namespace Tahoe 
#endif /* _HYPER_EVP_3D_ */
