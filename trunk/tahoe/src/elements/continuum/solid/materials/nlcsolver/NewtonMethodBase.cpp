/* $Id: NewtonMethodBase.cpp,v 1.5 2003-01-27 07:00:28 paklein Exp $ */
#include "NewtonMethodBase.h"
#include "MaterialsConfig.h"

#include "NLCSolver.h"
#include "dArrayT.h"
#include "LAdMatrixT.h"

#ifdef PLASTICITY_CRYSTAL_MATERIAL
#include "Utils.h"
#endif

using namespace Tahoe;

NewtonMethodBase::NewtonMethodBase() { }

void NewtonMethodBase::GetNewtonStep(NLCSolver& nlcsolve, dArrayT& step) const
{
  step.SetToScaled(-1., nlcsolve.fRHS);

  try { nlcsolve.fLHS.LinearSolve(step); }

  catch(ExceptionT::CodeT error)
    {
      if (NLCSolver::NLCS_MESSAGES) 
      {
#ifdef PLASTICITY_CRYSTAL_MATERIAL
         writeMessage("NewtonMethodBase::GetNewtonStep: Problems in LinearSolve");
#else
		ExceptionT::GeneralFail("NewtonMethodBase::GetNewtonStep", "Problems in LinearSolve");
#endif
	  }
      throw;
    }
}
