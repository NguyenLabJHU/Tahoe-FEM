/*
  File: EVPFDBaseT.h
*/

#ifndef _EVP_FD_BASE_T_H_
#define _EVP_FD_BASE_T_H_

#include "IsotropicT.h"
#include "FDHookeanMatT.h"

#include "NLCSolverWrapperPtr.h"
#include "KineticEqnBase.h"

#include "GlobalT.h"
#include "ArrayT.h"
#include "dArrayT.h"
#include "Array2DT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "ifstreamT.h" 

class NLCSolver;
class ElasticT;

class EVPFDBaseT : public FDHookeanMatT, public IsotropicT
{
 public:
  // constructor
  EVPFDBaseT(ifstreamT& in, const ElasticT& element);

  // destructor
  virtual ~EVPFDBaseT();

  // allocate space/initialize crystal arrays (all)
  virtual bool NeedsInitialization() const;
  virtual void Initialize();

  // required parameter flag
  virtual bool NeedLastDisp() const;

  // PVFs defined in derived classes 
  virtual void SetKineticEquation() = 0;
  virtual void InitializeVariables() = 0;
  virtual int  NumVariablesPerElement() = 0;

  // PVFs invoked from NLCSolver and defined in derived classes
  virtual void FormRHS(const dArrayT& variab, dArrayT& rhs) = 0;
  virtual void FormLHS(const dArrayT& variab, dMatrixT& lhs) = 0;

  // general accesors
  ifstreamT& Input();
  double Temperature();

  // print data read
  virtual void Print(ostream& out) const;

 protected:
  // print name
  virtual void PrintName(ostream& out) const;

  // size of system of nonlinear equations
  virtual int GetNumberOfEqns();

  // allocate all elements at once
  void AllocateElements();

  // deformation gradient
  virtual const dMatrixT& DeformationGradient(const LocalArrayT& disp);

 private:
  // solver for nonlinear constitutive equations
  void SetConstitutiveSolver();

 protected:
  // current time & time step
  const double& ftime;
  double fdt;

  // temperature
  double fTheta;

  // status flag
  const GlobalT::StateT& fStatus;

  // reference to displacements
  const LocalArrayT& fLocDisp;
  const LocalArrayT& fLocLastDisp;

  // stream for input data
  ifstreamT fInput;

  // Lame's elastic constants
  double fmu;
  double flambda;
  double fbulk;

  // control code for supporting classes
  int fSolverCode;
  int fKinEqnCode;

  // iteration counter of NLCSolver
  int fIterCount;

  // pointers to supporting classes
  KineticEqnBase* fKineticEqn;
  NLCSolver*  fSolver;

  // handle to NLCSolver
  NLCSolverWrapperPtr fSolverPtr;

  // total deformation gradient
  dMatrixT fFtot;

  // Cauchy stress
  dSymMatrixT fs_ij;

  // Tangent Moduli
  dMatrixT fc_ijkl;
};

inline ifstreamT& EVPFDBaseT::Input() { return fInput; }
inline double EVPFDBaseT::Temperature() { return fTheta; }

#endif /* _EVP_FD_BASE_T_H_ */
