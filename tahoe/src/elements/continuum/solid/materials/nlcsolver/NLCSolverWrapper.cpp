/*
  File: NLCSolverWrapper.cpp
*/

#include "NLCSolverWrapper.h"
#include "PolyCrystalMatT.h"
#include "SlipHardening.h"
#include "EVPFDBaseT.h"

/* base class: NLCSolverWrapper */

using namespace Tahoe;

NLCSolverWrapper::~NLCSolverWrapper() { }

/* derived class: SolverWrapperPoly */
SolverWrapperPoly::SolverWrapperPoly(PolyCrystalMatT& poly):
  fpoly(poly) { }

SolverWrapperPoly::~SolverWrapperPoly() { }

void SolverWrapperPoly::FormRHS(dArrayT& x, dArrayT& rhs) 
{ 
  fpoly.FormRHS(x, rhs); 
}

void SolverWrapperPoly::FormLHS(dArrayT& x, dMatrixT& lhs) 
{
  fpoly.FormLHS(x, lhs); 
}

/* derived class: SolverWrapperHard */
SolverWrapperHard::SolverWrapperHard(SlipHardening& hard):
  fhard(hard) { }

SolverWrapperHard::~SolverWrapperHard() { }

void SolverWrapperHard::FormRHS(dArrayT& x, dArrayT& rhs) 
{
  fhard.FormRHS(x, rhs); 
}

void SolverWrapperHard::FormLHS(dArrayT& x, dMatrixT& lhs)
{
  fhard.FormLHS(x, lhs); 
}

/* derived class: SolverWrapperEVPBase */
SolverWrapperEVPBase::SolverWrapperEVPBase(EVPFDBaseT& evp):
  fevp(evp) { }

SolverWrapperEVPBase::~SolverWrapperEVPBase() { }

void SolverWrapperEVPBase::FormRHS(dArrayT& x, dArrayT& rhs) 
{
  fevp.FormRHS(x, rhs); 
}

void SolverWrapperEVPBase::FormLHS(dArrayT& x, dMatrixT& lhs)
{
  fevp.FormLHS(x, lhs); 
}
