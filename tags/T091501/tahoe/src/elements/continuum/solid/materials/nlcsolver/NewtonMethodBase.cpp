/*
  File: NewtonMethodBase.cpp
*/

#include "NewtonMethodBase.h"
#include "NLCSolver.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"
#include "Utils.h"

NewtonMethodBase::NewtonMethodBase() { }

void NewtonMethodBase::GetNewtonStep(NLCSolver& nlcsolve, dArrayT& step) const
{
  step.SetToScaled(-1., nlcsolve.fRHS);

  try { nlcsolve.fLHS.LinearSolve(step); }

  catch(int error)
    {
      writeMessage("NewtonMethodBase::GetNewtonStep: Problems in LinearSolve");
      throw;
    }
}
