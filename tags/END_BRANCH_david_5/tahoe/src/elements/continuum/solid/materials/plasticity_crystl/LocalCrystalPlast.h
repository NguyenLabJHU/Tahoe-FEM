/* $Id: LocalCrystalPlast.h,v 1.10 2005-01-21 16:51:21 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_H_
#define _LOCAL_CRYSTAL_PLAST_H_

#include "PolyCrystalMatT.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "LAdMatrixT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;
class ElementCardT;
class StringT;

class LocalCrystalPlast : public PolyCrystalMatT
{
 public:
  // constructor
  LocalCrystalPlast(void);

  // destructor
  ~LocalCrystalPlast();

  // number of variables to compute in NLCSolver
  virtual int NumberOfUnknowns() const;

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

  // form of tangent matrix
  virtual GlobalT::SystemTypeT TangentType(void) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

 protected:

  // slip kinetics
  virtual void SetSlipKinetics();

  // slip hardening law
  virtual void SetSlipHardening();

  // recover crystal variables
  virtual void LoadCrystalData(ElementCardT& element, int intpt, int igrain);

  // recover averaged stress and moduli at IP
  virtual void LoadAggregateData(ElementCardT& element, int intpt);

  // solve for the state at each crystal
  virtual void IterateOnCrystalState(bool& stateConverged, int subIncr);

  // restores and saves solution for FDGamma (when subincrementation is on!)
  virtual void RestoreSavedSolution();
  virtual void SaveCurrentSolution();

  // trial deformation quantities
  void TrialDeformation();

  // forward gradient approximation for dgamma
  void ForwardGradientEstimate();

  // check convergence on DGamma
  bool Converged(double toler);

  // crystal Cauchy stress
  virtual void CrystalS_ij();

  // crystal consistent moduli
  virtual void CrystalC_ijkl();

 // elastic contribution to crystal moduli
  virtual void CrystalC_ijkl_Elastic();

  //plastic contribution to crystal moduli
  virtual void CrystalC_ijkl_Plastic();

  // solve for incremental shear strain
  virtual void SolveForDGamma() { };
  virtual void SolveForDGamma(int& ierr);

  // slip system resolve shear stress
  void ResolveShearStress();

  // inverse of incremental plastic deformation gradient
  virtual void DeltaFPInverse(const dArrayT& dgamma);

  // fDFpi/dDGamma term in local Jacobian
  void dDFpidDGamma(const dArrayT& dgamma, ArrayT<dMatrixT>& array);

  // frequent term used to compute the crystal moduli
  void dTaudCe(const dMatrixT& Z, const dSymMatrixT& P, dSymMatrixT& symmatx);

  // update plastic deformation gradient
  void FPInverse();

  // polar decomposition of deformation gradient
  void PolarDecomp();

 private:
 
   // number of crystal variables to be stored
  virtual int NumVariablesPerElement();

  // initial value of crystal variables
  virtual void InitializeCrystalVariables(ElementCardT& element);

 protected:
  // number of hardening variables (used in Grad derived class)
  int fNumHard;

  // norms to check convergence on DGamma
  double fMagDGam0;
  double fMagDGam;

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

  // crystal consistent tangent operator in Bbar
  dMatrixT fcBar_ijkl;  

  // Sym Schmidt tensors in sample coords
  ArrayT<dSymMatrixT> fP;

  // incremental slip system shearing rate
  dArrayT fDGamma_n;
  dArrayT fdgam_save;

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
  dSymMatrixT fsymmatx3;

  // workspaces for computing moduli 
  LAdMatrixT fLHS;
  ArrayT<dSymMatrixT> fA;
  ArrayT<dSymMatrixT> fB;
  ArrayT<dMatrixT> farray;

  // mid-point rule parameters
  double fTheta;

  // equivalent stress
  dSymMatrixT fAvgStress;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_H_ */
