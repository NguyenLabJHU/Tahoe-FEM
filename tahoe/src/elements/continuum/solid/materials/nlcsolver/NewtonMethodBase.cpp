/*
  File: NewtonMethodBase.cpp
*/

#include "NewtonMethodBase.h"
#include "NLCSolver.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"
#include "Utils.h"


using namespace Tahoe;

NewtonMethodBase::NewtonMethodBase() { }

void NewtonMethodBase::GetNewtonStep(NLCSolver& nlcsolve, dArrayT& step) const
{
  step.SetToScaled(-1., nlcsolve.fRHS);

  try { nlcsolve.fLHS.LinearSolve(step); }

  catch(ExceptionT::CodeT error)
    {
      if (NLCSolver::NLCS_MESSAGES) 
         writeMessage("NewtonMethodBase::GetNewtonStep: Problems in LinearSolve");
      throw;
    }
}
