/* $Id: NLCSolver.cpp,v 1.7 2003-01-29 07:35:00 paklein Exp $ */
#include "NLCSolver.h"

#include <iostream.h>
#include "NLCSolverWrapperPtr.h"

// some macros
#ifndef max
static int max(int i1, int i2) {return i1 >= i2 ? i1 : i2;};
static double max(double d1, double d2) {return d1 >= d2 ? d1 : d2;};
#endif

using namespace Tahoe;

const double STP_MAX = 1000.;
const double DBL_EPS = 2.2204460492503131E-16;

bool NLCSolver::NLCS_MESSAGES = false;

NLCSolver::NLCSolver(const int dim, const NewtonMethodBase& method)
  : // dimension of arrays
    fDim (dim),

    // default local Jacobian 
    fJacCode (kAnalyticalJ),
    //fJacCode (kFiniteDiffJ),
    //fJacCode (kSecantJ),

    // allocate space for arrays
    fGradF           (fDim),
    fNewtonStep      (fDim),
    fTypX            (fDim),
    fInvTypX         (fDim),
    fLastAcceptedX   (fDim),
    fLastAcceptedRHS (fDim),
    fRHS             (fDim),
    fLHS             (fDim),
    farray           (fDim),
  
    // default values
    fMaxConsecutiveMaxSteps (5),
    fMaxRejections          (10),
    fMaxIterations          (25),
    fFuncTol                (1.e-8),
    fDigits                 (-1),

    // solution method
    fMethod(method.clone())
{
  // default tolerance settings
  fStepTol = SetDefaultStepTol();
  fGradTol = sqrt(fStepTol);
}

NLCSolver::~NLCSolver() 
{
  delete fMethod;
}

void NLCSolver::Solve(NLCSolverWrapperPtr theModel, dArrayT& X, int& ierr)
{
  // reset some counters and bool variables
  ResetMemberData();

  // compute rhs and function value F
  theModel.FormRHS(X, fRHS);
  fF = 0.5 * dArrayT::Dot(fRHS, fRHS);
  
  // compute local Jacobian
  FormLocalJacobian(theModel, X);

  // set some default values
  SetDefaultTypValues(X);

  // check convergence of initial guess
  TestForInitialConvergence(X);

  if (NLCS_MESSAGES) {
    cout << "** NEWTON ITERATIONS **" << endl;
    cout << "   Iter = " << fIterationCount << endl;
  }

  // newton iterations
  while (!fConverged && !fFailed)
    {
      ComputeTrialPoint(X);

      theModel.FormRHS(X, fRHS);
      fF = 0.5 * dArrayT::Dot(fRHS, fRHS);

      TestTrialPoint(X);
	
      if (fTrialPointAccepted)
	{
	  FormLocalJacobian(theModel, X);

	  TestForConvergence(X);

	  if (NLCS_MESSAGES) {
	    cout << endl;
	    cout << "   iter = " << fIterationCount << endl; 
	  }
	}
    }

  // set returning code
  if (fConverged) ierr = 0;
  if (fFailed) ierr = 1;
}

void NLCSolver::SolveOneIteration(NLCSolverWrapperPtr theModel, dArrayT& X, int& ierr)
{
  // reset some counters and bool variables
  ResetMemberData();

  // compute rhs and function value F
  theModel.FormRHS(X, fRHS);
  fF = 0.5 * dArrayT::Dot(fRHS, fRHS);
  
  // check convergence of solution
  fConverged = FuncTest(X);
  if (fConverged)
    {
      ierr = 0;
      return;
    }
  fLastAcceptedX = X;
  fLastAcceptedF = fF;

  // compute local Jacobian
  FormLocalJacobian(theModel, X);

  // set some default values
  SetDefaultTypValues(X);

  // carry out one newton iteration
  while (!fFailed)
    {
      ComputeTrialPoint(X);

      theModel.FormRHS(X, fRHS);
      fF = 0.5 * dArrayT::Dot(fRHS, fRHS);

      TestTrialPoint(X);
	
      if (fTrialPointAccepted)
	{
	  ierr = 1;
	  return;
	}
    }

  // set code for failure
  ierr = 2;
  return;
}

/* called by external solve method */
void NLCSolver::SetRHS(const dArrayT& rhs)
{
  // assign rhs and compute function value F
  fRHS = rhs;
  fF = 0.5 * dArrayT::Dot(fRHS, fRHS);
}  

/* called by external solve method */
void NLCSolver::SetLHS(const dMatrixT& lhs)
{
  // assign lhs and compute gradient of F
  fLHS = lhs;
  fLHS.MultTx(fRHS, fGradF);
}

void NLCSolver::SetDefaultTypValues(const dArrayT& X)
{
  // set typical parameters and function value
  SetDefaultTypX(X);
  SetDefaultTypFunction(fF);

  // set max Newton's step
  SetDefaultMaxStep(X);
}

void NLCSolver::SetTypicalX(const dArrayT& typX)
{
	const char caller[] = "NLCSolver::SetTypicalX";

  if (typX.Length() != fDim) ExceptionT::GeneralFail(caller, "Dim(typX) != fDim");

  fTypX = typX;
  for (int i = 0; i < fDim; i++)
    {
      if (typX[i] <= 0.) ExceptionT::GeneralFail(caller, "typX_i < 0");
      fInvTypX[i] = 1./typX[i];
    }
}

void NLCSolver::ResetMemberData()
{
  fConverged = false;
  fFailed = false;
  fMaxStepTaken = false;
  fTrialPointAccepted = false;
  fRejectionCount = 0;
  fIterationCount = 1;
  fConsecutiveMaxSteps = 0;
}

void NLCSolver::TestForInitialConvergence(const dArrayT& X)
{
  // restricter test for initial guess
  fFuncTol *= 0.01;

  // check convergence on function value
  fConverged = FuncTest(X);
  if (fConverged && NLCS_MESSAGES)
    cout << "CONVERGENCE AT INITIAL GUESS." << endl;

  // reset function tolerance
  fFuncTol *= 100.;

  // if did not converge, save current solution and function value
  if (!fConverged) 
    {
      fLastAcceptedX = X;
      fLastAcceptedF = fF;
    }
}
	  
void NLCSolver::ComputeTrialPoint(dArrayT& X)  
{
#pragma unused(X)
  if (fRejectionCount == 0)
    {
      fMethod->GetNewtonStep(*this, fNewtonStep);
      farray = fNewtonStep;
      farray *= fInvTypX;
      fNewtonStepLength = farray.Magnitude();
      if (fNewtonStepLength > fMaxStep)
	{
	  fNewtonStep *= (fMaxStep / fNewtonStepLength);
	  fNewtonStepLength = fMaxStep;
	}
      if (NLCS_MESSAGES)
	cout << "Scaled Newton step length = " << fNewtonStepLength << endl;
    }
}

void NLCSolver::TestForConvergence(const dArrayT& X)
{
  if (FuncTest(X))
    {
      fConverged = true;
      if (NLCS_MESSAGES)
	cout << "SMALL DECREASE IN FUNCTION VALUE, ITERATION CONVERGED." << endl;
    }
  else if (StepTest(X))
    {
      fConverged = true;
      if (NLCS_MESSAGES)
	cout << "RELATIVE STEP CHANGE WITHIN TOLERANCE, POSSIBLE FAILURE." << endl;
    }
  else if (fConsecutiveMaxSteps >= fMaxConsecutiveMaxSteps)
    {
      fFailed = true;
      if (NLCS_MESSAGES)
      {
	cout << "TOOK MAXIMUM STEP " << fConsecutiveMaxSteps << " TIMES." << endl;
	cout << "FUNCTION MAY BE UNBOUNDED BELOW OR HAVE AN ASYMPTOTE." << endl;
      }
    }
  else if (fIterationCount > fMaxIterations)
    {
      fFailed = true;
      if (NLCS_MESSAGES)
	cout << "EXCEEDED MAXIMUM NUMBER OF ITERATIONS, ITERATION FAILED." << endl;
    }
  else if (GradTest(X))
    {
      fConverged = true;
      if (NLCS_MESSAGES)
	cout << "RELATIVE GRADIENT WITHIN TOLERANCE, POSSIBLE FAILURE." << endl;
    }
  else if (!fConverged)
    {
      fLastAcceptedX = X;
      fLastAcceptedF = fF;
    }
}

void NLCSolver::TestForFailure(const dArrayT& X)
{
#pragma unused(X)
  if (fRejectionCount > fMaxRejections)
    {
      fFailed = true;
      if (NLCS_MESSAGES)
	cout << "REJECTED TRIAL POINT " << fMaxRejections << " TIMES, ITERATION FAILED." << endl;
    }
}	

inline void NLCSolver::SetFDigits(const int digits)
{
  // set precision digits
  fDigits = digits;

  // reset defaults step and grad tolerances
  fStepTol = SetDefaultStepTol();
  fGradTol = sqrt(fStepTol);
}

void NLCSolver::Print(ostream& out) const
{
  // print solver data
  out << "       Max# iterations . . . . . . . . . . . . . = "
      << fMaxIterations << "\n";
  out << "       Tolerance convergence on F. . . . . . . . = "
      << fFuncTol << "\n";
  out << "       Tolerance convergence on GradF. . . . . . = "
      << fGradTol << "\n";
}

/* PROTECTED MEMBER FUNCTIONS */
void NLCSolver::AcceptTrialPoint(dArrayT& X)
{
#pragma unused(X)
  fIterationCount++;
  fRejectionCount = 0;
  fTrialPointAccepted = true;
  if (fMaxStepTaken) fConsecutiveMaxSteps++;
}

void NLCSolver::RejectTrialPoint(dArrayT& X)
{
  fRejectionCount++;
  fTrialPointAccepted = false;
  TestForFailure(X);
  X = fLastAcceptedX;
}

/* PRIVATE MEMBER FUNCTIONS */
double NLCSolver::SetDefaultStepTol()
{
  double eta;
  if (fDigits == -1)
    eta = DBL_EPS;
  else
    eta = max(DBL_EPS, pow(10.0, -fDigits));
  return eta;
}

void NLCSolver::SetDefaultTypX(const dArrayT& X)
{
  for (int i = 0; i < fDim; i++)
    {
      fTypX[i] = max(1., fabs(X[i]));
      fInvTypX[i] = 1. / fTypX[i];
    }
}

void NLCSolver::SetDefaultMaxStep(const dArrayT& X)
{
  farray = X;
  farray *= fInvTypX;
  fMaxStep = STP_MAX * max(farray.Magnitude(), fInvTypX.Magnitude());
}

void NLCSolver::SetDefaultTypFunction(const double& F)
{
  fTypF = max(1., fabs(F));
}

bool NLCSolver::FuncTest(const dArrayT& X) const
{
#pragma unused(X)
  double temp = 0.;
  for (int i = 0; i < fDim; i++)
    temp = max(fabs(fRHS[i]), temp);

  if (NLCS_MESSAGES) {
    cout << "NLCSolver::FuncTest: " << endl;
    cout << "   max(fRHS[i]) = " << temp << endl;
    cout << "   sqrt(2*fF)   = " << sqrt(2.*fF) << endl;
  }

  return (temp < fFuncTol);
}

bool NLCSolver::StepTest(const dArrayT& X) const 
{
  // test relative change in successive values of X
  double temp;
  int i = 0;
  bool test = true;
  while(i < fDim && test)
    {
      temp = (fabs(X[i] - fLastAcceptedX[i])) / max(fabs(X[i]), fTypX[i]);
      test = (temp < fStepTol);
      i++;
    }
  return test;
}

bool NLCSolver::GradTest(const dArrayT& X) const
{
  // test relative gradient
  double temp1 = max(fTypF, fabs(fF));
  int i = 0;
  bool test = true;
  double temp;
  while (i < fDim && test)
    {
      temp = fabs(fGradF[i]) * max(fabs(X[i]), fTypX[i]) / temp1;
      test = (temp < fGradTol);
      i++;
    }
  return test;
}

void NLCSolver::FormLocalJacobian(NLCSolverWrapperPtr theModel, dArrayT& X)
{
  // compute local Jacobian
  switch(fJacCode)
    {
    case kAnalyticalJ:
      theModel.FormLHS(X, fLHS);
      break;

    case kFiniteDiffJ:
      FiniteDifferenceJac(theModel, X);
      break;

    case kSecantJ:
      if (fIterationCount == 1) 
	FiniteDifferenceJac(theModel, X);   // initial jacobian
      else 
	SecantUpdateJac(theModel, X);
      break;

    default:
		ExceptionT::GeneralFail("NLCSolver::FormLocalJacobian", "Bad fJacCode: %d", fJacCode);
    }

  // compute gradient of F: GradF
  fLHS.MultTx(fRHS, fGradF);
}

void NLCSolver::FiniteDifferenceJac(NLCSolverWrapperPtr theModel, dArrayT& X)
{
  // relative noise
  // double eta = pow(10., -fDigits);
  double eta = DBL_EPS;
  eta = sqrt(eta);

  // zero out Jacobian
  fLHS = 0.;

  // forward difference approximation to local Jacobian
  for (int j = 0; j < fDim; j++)
    {
      // step size rule
      double stepsizej = eta * max(fabs(X[j]), fTypX[j]) * X[j]/fabs(X[j]);

      // to reduce finite precision errors slightly
      double tempj = X[j];
      X[j] += stepsizej;
      stepsizej = X[j] - tempj;

      // compute perturbed rhs
      theModel.FormRHS(X, farray);

      // column "j" of Jacobian
      for (int i = 0; i < fDim; i++)
	fLHS(i,j) = (farray[i] - fRHS[i]) / stepsizej;

      // reset Xj
      X[j] = tempj;
    }
}

void NLCSolver::SecantUpdateJac(NLCSolverWrapperPtr theModel, dArrayT& X)
{
  // relative noise
  // double eta = pow(10., -fDigits);
  double eta = DBL_EPS;

  // recompute previous rhs
  theModel.FormRHS(fLastAcceptedX, fLastAcceptedRHS);

  // difference between current and previous solution (s = X-Xlast)
  farray.SetToCombination(1., X, -1., fLastAcceptedX);

  // term: (s^T*Dx*Dx*s), Dx = 1/typX
  double denom = 0.;
  for (int i = 0; i < fDim; i++)
      denom += farray[i]*farray[i]*fInvTypX[i]*fInvTypX[i];

  // Broyden's update to Jacobian
  for (int i = 0; i < fDim; i++)
    {
      // term: (yi-Jik*sk) (use i-row)
      double tempi = fRHS[i] - fLastAcceptedRHS[i] - fLHS.DotRow(i, farray);

      // leave row-i unchanged if |(y-J*s)_i| is very small
      if (fabs(tempi) >= eta * (fabs(fRHS[i]) + fabs(fLastAcceptedRHS[i])))
	{
	  tempi /= denom;

	  // update i-row: Jij = Jij + (y-J*s)_i*(s*Dx*Dx)_j^T/(s^T*Dx*Dx*s)
	  for (int j = 0; j < fDim; j++)
	    fLHS(i,j) += tempi * farray[j] * (fInvTypX[j]*fInvTypX[j]);
	}
    }
}
