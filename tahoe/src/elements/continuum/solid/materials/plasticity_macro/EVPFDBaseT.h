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

	/* set (material) tangent modulus */
	virtual void SetModulus(dMatrixT& modulus);

  // print name
  virtual void PrintName(ostream& out) const;

  // size of system of nonlinear equations
  virtual int GetNumberOfEqns();

  // allocate all elements at once
  void AllocateElements();

  // deformation gradient
//  virtual const dMatrixT& DeformationGradient(const LocalArrayT& disp);
//DEV - let's see who needs this

	// function to compute 3D deformation regardless of dimensionality of the
	// problem. For 2D, the out-of-plane direction is x3 and the deformation
	// is assumed to be plane strain
	void Compute_Ftot_3D(dMatrixT& F_3D) const;	

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

  // 3D total deformation gradient
  dMatrixT fFtot;

  // Cauchy stress
  dSymMatrixT fs_ij;

  // Tangent Moduli
  dMatrixT fc_ijkl;
};

inline ifstreamT& EVPFDBaseT::Input() { return fInput; }
inline double EVPFDBaseT::Temperature() { return fTheta; }

#endif /* _EVP_FD_BASE_T_H_ */
