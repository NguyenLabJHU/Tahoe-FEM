/* $Id: LocalCrystalPlastFp.h,v 1.6 2004-07-15 08:29:07 paklein Exp $ */
#ifndef _LOCAL_CRYSTAL_PLAST_FP_H_
#define _LOCAL_CRYSTAL_PLAST_FP_H_

#include "PolyCrystalMatT.h"

#include <iostream.h>
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "LAdMatrixT.h"
#include "SpectralDecompT.h"

namespace Tahoe {

class ifstreamT;
class SolidElementT;
class ElementCardT;
class StringT;

class LocalCrystalPlastFp : public PolyCrystalMatT
{
 public:
  // constructor
  LocalCrystalPlastFp(ifstreamT& in, const FSMatSupportT& support);

  // destructor
  ~LocalCrystalPlastFp();

  // number of variables to compute in NLCSolver
  virtual int NumberOfUnknowns() const;

  // Cauchy stress - Taylor average    
  virtual const dSymMatrixT& s_ij();   

  // modulus - Taylor average 
  virtual const dMatrixT& c_ijkl();

  // form residual
  virtual void FormRHS(const dArrayT& fparray, dArrayT& rhs);

  // form Jacobian
  virtual void FormLHS(const dArrayT& fparray, dMatrixT& lhs);

  // update/reset crystal state
  virtual void UpdateHistory();
  virtual void ResetHistory();

  // output related methods
  virtual int NumOutputVariables() const;
  virtual void OutputLabels(ArrayT<StringT>& labels) const;
  virtual void ComputeOutput(dArrayT& output);

  // form of tangent matrix
  virtual GlobalT::SystemTypeT TangentType(void) const;

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

  // solution related methods for Fp and hardness
  virtual void InitialEstimateForFp();
  virtual void InitialEstimateForHardening();
  virtual void SolveForPlasticDefGradient(int& ierr);
  virtual void SolveForHardening();

  // restores and saves solution for Fp (used when subincrementation is on!)
  virtual void RestoreSavedSolution();
  virtual void SaveCurrentSolution();

  // check convergence on Fp
  bool Converged(double toler);

  // crystal Cauchy stress
  virtual void CrystalS_ij();

  // crystal consistent moduli
  virtual void CrystalC_ijkl();

  // compute modified Schmidt tensor : X = fLHS^(-1) * Z
  void ComputeSchmidtXTensor();

  // slip system resolve shear stress
  virtual void ResolveShearStress();

  // switches between 3x3 matrix (tensor) and 9x1 array representations: A <-> { A }
  void Rank2ToArray9x1(const dMatrixT& matrix, dArrayT& array);
  void Rank2FromArray9x1(dMatrixT& matrix, const dArrayT& array);

  // outer product of 9x1 arrays
  void OuterProduct9x9(const dMatrixT& matrix1, const dMatrixT& matrix2, dMatrixT& outer);

  // constructs 9x9 matrices from 3x3 matrices: X*A=[AR]{X}, A*X=[AL]{X}
  void Rank2ToMatrixAR9x9(const dMatrixT& matrix, dMatrixT& AR);
  void Rank2ToMatrixAL9x9(const dMatrixT& matrix, dMatrixT& AL);

  // add gradient dependent term for Local Jacobian 
  virtual void AddGradTermToLHS(dMatrixT& lhs, const dMatrixT& matx);

  // add gradient dependent term to consistent moduli
  virtual void AddGradTermToC_ijkl();

 private:
 
   // number of crystal variables to be stored
  virtual int NumVariablesPerElement(void);

  // initial value of crystal variables
  virtual void InitializeCrystalVariables(ElementCardT& element);

 protected:
  // penalty parameter for detFp
  double fPenalty;

  // norms to check convergence on Fp during state iters
  double fFpNorm0;
  double fFpNorm;

  // elastic deformation gradients
  dMatrixT fFe_n;
  dMatrixT fFe;

  // plastic deformation gradients
  dMatrixT fFp_n;
  dMatrixT fFp;
  dMatrixT fFpi;
  dMatrixT fFp_save;

  // right Cauchy-Green tensor
  dSymMatrixT fC;

  // tensors in intermediate configuration
  dSymMatrixT fCeBar;
  dSymMatrixT fEeBar;
  dSymMatrixT fSBar;

  // crystal consistent tangent operator (in Bbar)
  dMatrixT fcBar_ijkl;  

  // work spaces used in polar decomposition of fFe
  dArrayT fEigs;
  dMatrixT fRe;
  dSymMatrixT fUe;
  SpectralDecompT fSpecD;

  // second order identity tensor
  dMatrixT fIMatx;
  dSymMatrixT fISym;

  // workspaces: 3x3 (sym & unsym) and 6x6 matrices
  dSymMatrixT fSymMatx1;
  dSymMatrixT fSymMatx2;
  dSymMatrixT fSymMatx3;
  dMatrixT fMatxCe;
  dMatrixT fMatxSb;
  dMatrixT fMatx1;
  dMatrixT fMatx2;
  dMatrixT fMatx3;
  dMatrixT fRank4;

  // work spaces: 9x1 arrays,
  dArrayT fFpArray;
  dArrayT fArray1;
  dArrayT fArray2;

  // work spaces: 9x9 matrices
  dMatrixT fRankIV_1;
  dMatrixT fRankIV_2;
  dMatrixT fRankIV_3;
  dMatrixT fLHS;

  // work spaces: arrays of 3x3 matrices
  ArrayT<dSymMatrixT> fA;
  ArrayT<dSymMatrixT> fB;
  ArrayT<dMatrixT> fArrayOfMatx;

  // equivalent stress
  dSymMatrixT fAvgStress;
};

} // namespace Tahoe 
#endif /* _LOCAL_CRYSTAL_PLAST_FP_H_ */

