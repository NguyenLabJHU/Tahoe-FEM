/*
  File: NLCSolver.h
*/

#ifndef _NLC_SOLVER_H_
#define _NLC_SOLVER_H_

#include <iostream.h>

#include "StringT.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"
#include "NewtonMethodBase.h"


namespace Tahoe {

class NLCSolverWrapperPtr;

class NLCSolver
{
  friend class NewtonMethodBase;

 public:
  // enum variable for solver type
  enum SolverType { kNLCSolver_LS = 1,
		    kNLCSolver_TR = 2 };

  // enum variable for local Jacobian computation
  enum JacType { kAnalyticalJ = 1,
		 kFiniteDiffJ = 2,
		 kSecantJ     = 3 };

  // constructor
  NLCSolver(const int dim, const NewtonMethodBase& method);

  // destructor
  virtual ~NLCSolver();

  // embedded solve methods for nonlinear equations (Newton's method)
  virtual void Solve(NLCSolverWrapperPtr theModel, dArrayT& X, int& ierr);
  virtual void SolveOneIteration(NLCSolverWrapperPtr theModel, dArrayT& X, int& ierr);

  // set rhs and lhs of system Ax = b (lhs=A, rhs=b)
  void SetRHS(const dArrayT& rhs);
  void SetLHS(const dMatrixT& lhs);

  // set typical (scaling) values and max Newton's step
  void SetDefaultTypValues(const dArrayT& X);

  // accessors to typical parameters (X) values
  void SetTypicalX(const dArrayT& typX);
  dArrayT& GetTypicalX();

  // reinitialize counters and bool variables
  virtual void ResetMemberData();

  // test convergence of intial guess
  void TestForInitialConvergence(const dArrayT& X);

  // compute/test trial solution
  virtual void ComputeTrialPoint(dArrayT& X);
  virtual void TestTrialPoint(dArrayT& X) = 0;

  // tests for convergence
  void TestForConvergence(const dArrayT& X);
  void TestForFailure(const dArrayT& X);

  // accesors to set/get stopping tolerances
  void SetStepTol(double tol);
  void SetFuncTol(double tol);
  void SetGradTol(double tol);
  double GetStepTol() const;
  double GetFuncTol() const;
  double GetGradTol() const;

  // accesors to set/get max values of step and counters
  void SetMaxStep(double step);
  void SetMaxConsecutiveMaxSteps(int n);   // default 5
  void SetMaxRejections(int n);            // default 10
  void SetMaxIterations(int n);            // default 25
  double GetMaxStep() const;
  int GetMaxConsecutiveMaxSteps() const;
  int GetMaxRejections() const;
  int GetMaxIterations() const;

  // accesors to counters
  int GetIterationCount() const;
  int GetRejectionCount() const;

  // accesor to set the solution method
  void SetMethod(const NewtonMethodBase& method);

  // additional accesors
  bool TrialPointAccepted() const;
  bool Converged() const;
  bool Failed() const;

  // access to newton step
  void GetNewtonStep(dArrayT& step) const;

  // accessors to set/get precision digits
  void SetFDigits(const int digits);      //will also reset step and grad tol's
  int GetFDigits() const;

  // set local Jacobian code
  void SetJacobianCode(const int code);

  // print solver data
  virtual void Print(ostream& out) const;

 protected:
  // computations for trial solution
  void AcceptTrialPoint(dArrayT& X);
  void RejectTrialPoint(dArrayT& X);

 private:
  // set some default values
  double SetDefaultStepTol();
  void SetDefaultTypX(const dArrayT& X);
  void SetDefaultMaxStep(const dArrayT& X);
  void SetDefaultTypFunction(const double& F);

  // checks of step, function and gradient values
  bool FuncTest(const dArrayT& X) const;
  bool GradTest(const dArrayT& X) const;
  bool StepTest(const dArrayT& X) const;

  // local Jacobian computation
  void FormLocalJacobian(NLCSolverWrapperPtr theModel, dArrayT& X);
  void FiniteDifferenceJac(NLCSolverWrapperPtr theModel, dArrayT& X);
  void SecantUpdateJac(NLCSolverWrapperPtr theModel, dArrayT& X);

 public:
  // bool variable to print messages
  static bool NLCS_MESSAGES;

 protected:
  // trial-point related bool variable 
  bool fTrialPointAccepted;

  // dimension of the problem
  int fDim;

  // local Jacobian computation code
  int fJacCode;

  // number of significative digits in solution
  int fDigits;

  // counters
  int fRejectionCount;
  int fIterationCount;

  // function-related variables
  double fF;
  double fTypF;
  double fLastAcceptedF;
  dArrayT fGradF;

  // step-related variables
  bool fMaxStepTaken;
  double fMaxStep;
  double fNewtonStepLength;
  double fStepTol;
  dArrayT fNewtonStep;  

  // typical values for unknown (scaling)
  dArrayT fTypX;
  dArrayT fInvTypX;
  
  // last value of unknown
  dArrayT fLastAcceptedX;
  dArrayT fLastAcceptedRHS;

  // lhs and rhs of system: Ax=b (lhs=A, rhs=b)
  dArrayT fRHS;
  LAdMatrixT fLHS;

  // workspace
  dArrayT farray;

  // solution method
  NewtonMethodBase* fMethod;

 private:
  // bool variables to check convergence/failure status of solution
  bool fConverged;
  bool fFailed;

  // counter for max step taken
  int fConsecutiveMaxSteps;

  // max values of counters
  int fMaxConsecutiveMaxSteps;
  int fMaxRejections;
  int fMaxIterations;

  // tolerance for convergence checks
  double fFuncTol;
  double fGradTol;
};

inline dArrayT& NLCSolver::GetTypicalX() { return fTypX; }

inline void NLCSolver::SetStepTol(double tol) { fStepTol = tol; }
inline void NLCSolver::SetFuncTol(double tol) { fFuncTol = tol; }
inline void NLCSolver::SetGradTol(double tol) { fGradTol = tol; }
inline double NLCSolver::GetStepTol() const { return fStepTol; }
inline double NLCSolver::GetFuncTol() const { return fFuncTol; }
inline double NLCSolver::GetGradTol() const { return fGradTol; }

inline void NLCSolver::SetMaxStep(double step) { fMaxStep = step; }
inline void NLCSolver::SetMaxConsecutiveMaxSteps(int n) { fMaxConsecutiveMaxSteps = n; }
inline void NLCSolver::SetMaxRejections(int n) { fMaxRejections = n; }
inline void NLCSolver::SetMaxIterations(int n) { fMaxIterations = n; }

inline double NLCSolver::GetMaxStep() const { return fMaxStep; }
inline int NLCSolver::GetMaxConsecutiveMaxSteps() const { return fMaxConsecutiveMaxSteps; }
inline int NLCSolver::GetMaxRejections() const { return fMaxRejections; }
inline int NLCSolver::GetMaxIterations() const { return fMaxIterations; }

inline int NLCSolver::GetIterationCount() const { return fIterationCount; }
inline int NLCSolver::GetRejectionCount() const { return fRejectionCount; }

inline void NLCSolver::SetMethod(const NewtonMethodBase& method)
{
  if (fMethod) delete fMethod;
  fMethod = method.clone();
};

inline bool NLCSolver::TrialPointAccepted() const { return fTrialPointAccepted; }
inline bool NLCSolver::Converged() const { return fConverged; }
inline bool NLCSolver::Failed() const { return fFailed; }
inline void NLCSolver::GetNewtonStep(dArrayT& step) const { step = fNewtonStep; }

inline int NLCSolver::GetFDigits() const { return fDigits; }

inline void NLCSolver::SetJacobianCode(const int code) { fJacCode = code; }

} // namespace Tahoe 
#endif /*_NLC_SOLVER_H_ */





